use super::*;
use crate::logger::{append_json_log_fields, append_log_fields};
use crate::ordered_pipeline::{IndexedResult, PipelineConfig, run_bounded_ordered_pipeline};

pub(super) struct Processor {
    pub matcher: PatternMatcher,
    pub patterns: Vec<String>,
    pub filename: String,
    pub tag: [u8; 2],
    pub plain_log: bool,
    pub json_log: bool,
    pub filter: bool,
    pub invert: bool,
    pub suppress: bool,
}

pub(super) struct Summary {
    pub records: usize,
    pub bases: usize,
    pub hits: usize,
    pub records_hit: usize,
    pub counts: Vec<u32>,
}

impl Summary {
    pub fn new(patterns: usize) -> Self {
        Self {
            records: 0,
            bases: 0,
            hits: 0,
            records_hit: 0,
            counts: vec![0; patterns],
        }
    }
}

struct Work {
    index: u64,
    records: Vec<bam::Record>,
}

struct Output {
    records: Vec<bam::Record>,
    plain: String,
    json: Vec<u8>,
    summary: Summary,
}

// Keep ordinary processing errors in the ordered stream so the unchanged
// pipeline can drain and join all threads before we return the error.
struct BatchResult {
    index: u64,
    count: u64,
    result: Result<Output>,
}

impl IndexedResult for BatchResult {
    fn index(&self) -> u64 {
        self.index
    }
    fn index_span(&self) -> u64 {
        self.count
    }
}

fn read_chunk(reader: &mut dyn RecordReader, index: u64, size: usize) -> Result<Option<Work>> {
    let mut records = Vec::with_capacity(size);
    for _ in 0..size {
        let mut record = bam::Record::new();
        if !reader
            .read_into(&mut record)
            .context("Error reading SAM/BAM record")?
        {
            break;
        }
        records.push(record);
    }
    Ok((!records.is_empty()).then_some(Work { index, records }))
}

impl Processor {
    fn process(&self, work: Work) -> Result<Output> {
        let logging = self.plain_log || self.json_log;
        let mut output = Output {
            records: Vec::new(),
            plain: String::new(),
            json: Vec::new(),
            summary: Summary::new(if logging { self.patterns.len() } else { 0 }),
        };
        let mut json_first = true;
        let mut matched = vec![false; self.patterns.len()];
        for mut record in work.records {
            matched.fill(false);
            let sequence = record.sequence().to_vec();
            if logging {
                output.summary.records += 1;
                output.summary.bases += record.query_len() as usize;
                self.matcher.for_each_match(&sequence, |hit| {
                    matched[hit.pattern_index] = true;
                    output.summary.hits += 1;
                    output.summary.counts[hit.pattern_index] += 1;
                    let pattern = &self.patterns[hit.pattern_index];
                    if self.plain_log {
                        append_log_fields(
                            &mut output.plain,
                            &self.filename,
                            record.name(),
                            pattern,
                            hit.position,
                        );
                    }
                    if self.json_log {
                        append_json_log_fields(
                            &mut output.json,
                            &mut json_first,
                            &self.filename,
                            record.name(),
                            pattern,
                            hit.position,
                        );
                    }
                });
            } else {
                self.matcher.mark_matching_patterns(&sequence, &mut matched);
            }
            let has_match = matched.iter().any(|&value| value);
            if logging && has_match {
                output.summary.records_hit += 1;
            }
            if (self.filter && !has_match) || (self.invert && has_match) {
                continue;
            }
            let mut kmers: Vec<String> = self
                .patterns
                .iter()
                .zip(&matched)
                .filter(|(_, found)| **found)
                .map(|(pattern, _)| pattern.clone())
                .collect();
            match record.tags().get(&self.tag) {
                Some(tags::TagValue::String([], _)) | None => (),
                Some(tags::TagValue::String(value, _)) => {
                    kmers.extend(
                        from_utf8(value)
                            .context("Error reading existing tag value as UTF-8")?
                            .split(',')
                            .map(String::from),
                    );
                }
                _ => anyhow::bail!("Invalid tag value format. Expected string value."),
            }
            kmers.sort_unstable();
            kmers.dedup();
            record
                .tags_mut()
                .push_string(&self.tag, kmers.join(",").as_bytes());
            if !self.suppress {
                output.records.push(record);
            }
        }
        Ok(output)
    }
}

pub(super) fn run(
    mut reader: Box<dyn RecordReader + Send>,
    processor: Processor,
    threads: usize,
    chunk_size: usize,
    mut writer: Option<Box<dyn RecordWriter>>,
    logger: &mut BufferedLogger,
    json_logger: &mut Option<JsonLogger>,
) -> Result<Summary> {
    let mut summary = Summary::new(if processor.plain_log || processor.json_log {
        processor.patterns.len()
    } else {
        0
    });
    let mut consume = |output: Output| -> Result<()> {
        if let Some(writer) = writer.as_mut() {
            for record in &output.records {
                writer
                    .write(record)
                    .context("Error writing tagged record")?;
            }
        }
        logger.log_fragment(&output.plain)?;
        if let Some(logger) = json_logger.as_mut() {
            logger.log_fragment(&output.json)?;
        }
        summary.records += output.summary.records;
        summary.bases += output.summary.bases;
        summary.hits += output.summary.hits;
        summary.records_hit += output.summary.records_hit;
        for (total, count) in summary.counts.iter_mut().zip(output.summary.counts) {
            *total += count;
        }
        Ok(())
    };
    match threads {
        1 => {
            let mut index = 0;
            while let Some(work) = read_chunk(reader.as_mut(), index, chunk_size)? {
                index += work.records.len() as u64;
                consume(processor.process(work)?)?;
            }
        }
        2 => {
            // Combine reading and matching on one background thread; the main
            // thread writes. The generic pipeline requires at least three threads.
            let (tx, rx) = crossbeam_channel::bounded(2);
            let producer = std::thread::spawn(move || -> Result<()> {
                let mut index = 0;
                while let Some(work) = read_chunk(reader.as_mut(), index, chunk_size)? {
                    index += work.records.len() as u64;
                    if tx.send(processor.process(work)?).is_err() {
                        break;
                    }
                }
                Ok(())
            });
            let consumed = rx.iter().try_for_each(&mut consume);
            drop(rx);
            let produced = producer
                .join()
                .map_err(|_| anyhow::anyhow!("Tag processing thread panicked"))?;
            consumed?;
            produced?;
        }
        _ => {
            // Budget includes the producer and the main-thread consumer.
            let mut first_error = None;
            let pipeline_result = run_bounded_ordered_pipeline(
                PipelineConfig::new(threads - 2),
                move |tx| {
                    let mut index = 0;
                    while let Some(work) = read_chunk(reader.as_mut(), index, chunk_size)? {
                        index += work.records.len() as u64;
                        tx.send(work)
                            .map_err(|_| anyhow::anyhow!("Tag work queue closed"))?;
                    }
                    Ok(())
                },
                move |work: Work| {
                    Ok(BatchResult {
                        index: work.index,
                        count: work.records.len() as u64,
                        result: processor.process(work),
                    })
                },
                |batch| {
                    if first_error.is_none() {
                        if let Err(error) = batch.result.and_then(&mut consume) {
                            first_error = Some(error);
                        }
                    }
                    Ok(())
                },
            );
            if let Some(error) = first_error {
                return Err(error);
            }
            pipeline_result?;
        }
    }
    if let Some(writer) = writer.as_mut() {
        writer.finish().context("Error finishing tagged output")?;
    }
    Ok(summary)
}
