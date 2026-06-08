use crate::fastx_output::FastxFormat;
use crate::pattern_matching::PatternMatcher;
use crossbeam_channel::bounded;
use std::collections::BTreeMap;
use std::sync::Arc;
use std::thread;

#[derive(Debug, Clone, Copy, Eq, PartialEq)]
pub enum FileSlot {
    SingleOrFirst,
    Second,
}

#[derive(Debug, Clone, Eq, PartialEq)]
pub struct OutputRecord {
    pub id: Vec<u8>,
    pub seq: Vec<u8>,
    pub qual: Option<Vec<u8>>,
    pub format: FastxFormat,
}

#[derive(Debug, Clone, Eq, PartialEq)]
pub struct RecordHit {
    pub file_slot: FileSlot,
    pub record_id: Vec<u8>,
    pub pattern_index: usize,
    pub position: usize,
}

#[derive(Debug, Clone, Eq, PartialEq)]
pub struct RecordStats {
    pub num_bases: usize,
    pub hit_count: usize,
    pub distinct_record_hit: bool,
    pub pattern_hit_counts: Vec<u32>,
}

impl RecordStats {
    fn new(num_bases: usize, pattern_count: usize) -> Self {
        Self {
            num_bases,
            hit_count: 0,
            distinct_record_hit: false,
            pattern_hit_counts: vec![0; pattern_count],
        }
    }
}

#[derive(Debug, Clone, Eq, PartialEq)]
pub struct SingleResult {
    pub record_index: u64,
    pub matched: bool,
    pub extracted: bool,
    pub output: Option<OutputRecord>,
    pub stats: RecordStats,
    pub hits: Vec<RecordHit>,
}

#[derive(Debug, Clone, Eq, PartialEq)]
pub struct PairedResult {
    pub pair_index: u64,
    pub matched: bool,
    pub extracted: bool,
    pub output_1: Option<OutputRecord>,
    pub output_2: Option<OutputRecord>,
    pub stats_1: RecordStats,
    pub stats_2: RecordStats,
    pub hits: Vec<RecordHit>,
}

#[derive(Debug, Clone, Copy)]
pub struct RecordInput<'a> {
    pub id: &'a [u8],
    pub seq: &'a [u8],
    pub qual: Option<&'a [u8]>,
    pub format: FastxFormat,
}

#[derive(Debug, Clone)]
pub struct ExtractProcessor {
    pub matcher: Arc<PatternMatcher>,
    pub pattern_count: usize,
    pub logging_active: bool,
    pub suppress_output: bool,
    pub invert_match: bool,
}

impl ExtractProcessor {
    pub fn new(
        matcher: Arc<PatternMatcher>,
        pattern_count: usize,
        logging_active: bool,
        suppress_output: bool,
        invert_match: bool,
    ) -> Self {
        Self {
            matcher,
            pattern_count,
            logging_active,
            suppress_output,
            invert_match,
        }
    }

    pub fn process_single_record(
        &self,
        record_index: u64,
        record: RecordInput<'_>,
    ) -> SingleResult {
        let mut stats = RecordStats::new(record.seq.len(), self.pattern_count);
        let mut hits = Vec::new();
        let matched =
            self.process_record_matches(FileSlot::SingleOrFirst, record, &mut stats, &mut hits);
        let extracted = matched != self.invert_match;
        let output = self.output_record_if_needed(extracted, record);

        SingleResult {
            record_index,
            matched,
            extracted,
            output,
            stats,
            hits,
        }
    }

    pub fn process_paired_record(
        &self,
        pair_index: u64,
        read_1: RecordInput<'_>,
        read_2: RecordInput<'_>,
    ) -> PairedResult {
        let mut stats_1 = RecordStats::new(read_1.seq.len(), self.pattern_count);
        let mut stats_2 = RecordStats::new(read_2.seq.len(), self.pattern_count);
        let mut hits = Vec::new();

        let matched_1 =
            self.process_record_matches(FileSlot::SingleOrFirst, read_1, &mut stats_1, &mut hits);
        let matched_2 = if !self.logging_active && matched_1 {
            false
        } else {
            self.process_record_matches(FileSlot::Second, read_2, &mut stats_2, &mut hits)
        };
        let matched = matched_1 || matched_2;
        let extracted = matched != self.invert_match;
        let output_1 = self.output_record_if_needed(extracted, read_1);
        let output_2 = self.output_record_if_needed(extracted, read_2);

        PairedResult {
            pair_index,
            matched,
            extracted,
            output_1,
            output_2,
            stats_1,
            stats_2,
            hits,
        }
    }

    fn process_record_matches(
        &self,
        file_slot: FileSlot,
        record: RecordInput<'_>,
        stats: &mut RecordStats,
        hits: &mut Vec<RecordHit>,
    ) -> bool {
        if self.logging_active {
            self.matcher.for_each_match(record.seq, |hit| {
                stats.hit_count += 1;
                stats.distinct_record_hit = true;
                stats.pattern_hit_counts[hit.pattern_index] += 1;
                hits.push(RecordHit {
                    file_slot,
                    record_id: record.id.to_vec(),
                    pattern_index: hit.pattern_index,
                    position: hit.position,
                });
            });
            stats.distinct_record_hit
        } else {
            self.matcher.find_any(record.seq)
        }
    }

    fn output_record_if_needed(
        &self,
        extracted: bool,
        record: RecordInput<'_>,
    ) -> Option<OutputRecord> {
        if self.suppress_output || !extracted {
            return None;
        }

        Some(OutputRecord {
            id: record.id.to_vec(),
            seq: record.seq.to_vec(),
            qual: record.qual.map(|qual| qual.to_vec()),
            format: record.format,
        })
    }
}

pub trait IndexedResult {
    fn index(&self) -> u64;
}

impl IndexedResult for SingleResult {
    fn index(&self) -> u64 {
        self.record_index
    }
}

impl IndexedResult for PairedResult {
    fn index(&self) -> u64 {
        self.pair_index
    }
}

#[derive(Debug)]
pub struct OrderedResultBuffer<T> {
    next_expected_index: u64,
    pending: BTreeMap<u64, T>,
}

impl<T> OrderedResultBuffer<T> {
    pub fn new() -> Self {
        Self::with_start_index(0)
    }

    pub fn with_start_index(next_expected_index: u64) -> Self {
        Self {
            next_expected_index,
            pending: BTreeMap::new(),
        }
    }

    pub fn next_expected_index(&self) -> u64 {
        self.next_expected_index
    }

    pub fn pending_len(&self) -> usize {
        self.pending.len()
    }
}

impl<T: IndexedResult> OrderedResultBuffer<T> {
    pub fn push(&mut self, result: T) -> Option<T> {
        self.pending.insert(result.index(), result)
    }

    pub fn drain_ready(&mut self) -> Vec<T> {
        let mut ready = Vec::new();
        while let Some(result) = self.pending.remove(&self.next_expected_index) {
            self.next_expected_index += 1;
            ready.push(result);
        }
        ready
    }
}

impl<T> Default for OrderedResultBuffer<T> {
    fn default() -> Self {
        Self::new()
    }
}

#[derive(Debug, Clone, Eq, PartialEq)]
pub struct ExtractSummary {
    pub nb_records_tot: usize,
    pub nb_bases: usize,
    pub nb_hits_tot: (usize, usize),
    pub nb_records_hit: (usize, usize),
    pub nb_records_extracted: usize,
    pub pattern_hit_counts: Vec<u32>,
}

#[derive(Debug, Clone, Copy, Eq, PartialEq)]
pub struct PipelineConfig {
    pub worker_count: usize,
    pub work_queue_bound: usize,
    pub result_queue_bound: usize,
}

impl PipelineConfig {
    pub fn new(worker_count: usize) -> Self {
        Self {
            worker_count,
            work_queue_bound: worker_count.saturating_mul(2).max(1),
            result_queue_bound: worker_count.saturating_mul(2).max(1),
        }
    }
}

pub fn run_bounded_ordered_pipeline<W, R, P, C>(
    work_items: Vec<W>,
    config: PipelineConfig,
    process_work: P,
    mut consume_ready: C,
) -> anyhow::Result<()>
where
    W: Send + 'static,
    R: IndexedResult + Send + 'static,
    P: Fn(W) -> anyhow::Result<Vec<R>> + Send + Sync + 'static,
    C: FnMut(R) -> anyhow::Result<()>,
{
    if config.worker_count == 0 {
        anyhow::bail!("Pipeline worker count must be greater than zero.");
    }
    if config.work_queue_bound == 0 {
        anyhow::bail!("Pipeline work queue bound must be greater than zero.");
    }
    if config.result_queue_bound == 0 {
        anyhow::bail!("Pipeline result queue bound must be greater than zero.");
    }

    let (work_tx, work_rx) = bounded::<W>(config.work_queue_bound);
    let (result_tx, result_rx) = bounded::<R>(config.result_queue_bound);

    let producer = thread::spawn(move || -> anyhow::Result<()> {
        for work in work_items {
            work_tx.send(work).map_err(|_| {
                anyhow::anyhow!("Pipeline work queue closed before all work was sent.")
            })?;
        }
        Ok(())
    });

    let process_work = Arc::new(process_work);
    let mut workers = Vec::with_capacity(config.worker_count);
    for _ in 0..config.worker_count {
        let work_rx = work_rx.clone();
        let result_tx = result_tx.clone();
        let process_work = Arc::clone(&process_work);
        workers.push(thread::spawn(move || -> anyhow::Result<()> {
            while let Ok(work) = work_rx.recv() {
                for result in process_work(work)? {
                    result_tx.send(result).map_err(|_| {
                        anyhow::anyhow!(
                            "Pipeline result queue closed before all results were sent."
                        )
                    })?;
                }
            }
            Ok(())
        }));
    }
    drop(result_tx);

    let mut ordered_results = OrderedResultBuffer::new();
    let mut consume_error = None;
    while let Ok(result) = result_rx.recv() {
        ordered_results.push(result);
        for ready in ordered_results.drain_ready() {
            if let Err(error) = consume_ready(ready) {
                consume_error = Some(error);
                break;
            }
        }
        if consume_error.is_some() {
            break;
        }
    }
    drop(result_rx);

    let producer_result = producer
        .join()
        .map_err(|_| anyhow::anyhow!("Pipeline producer thread panicked."))?;

    let mut worker_error = None;
    for worker in workers {
        match worker.join() {
            Ok(Ok(())) => {}
            Ok(Err(error)) if worker_error.is_none() => worker_error = Some(error),
            Ok(Err(_)) => {}
            Err(_) if worker_error.is_none() => {
                worker_error = Some(anyhow::anyhow!("Pipeline worker thread panicked."));
            }
            Err(_) => {}
        }
    }

    producer_result?;
    if let Some(error) = worker_error {
        return Err(error);
    }
    if let Some(error) = consume_error {
        return Err(error);
    }
    if ordered_results.pending_len() > 0 {
        anyhow::bail!(
            "Pipeline completed with {} out-of-order result(s) still pending; missing result index {}.",
            ordered_results.pending_len(),
            ordered_results.next_expected_index()
        );
    }

    Ok(())
}

impl ExtractSummary {
    pub fn new(pattern_count: usize) -> Self {
        Self {
            nb_records_tot: 0,
            nb_bases: 0,
            nb_hits_tot: (0, 0),
            nb_records_hit: (0, 0),
            nb_records_extracted: 0,
            pattern_hit_counts: vec![0; pattern_count],
        }
    }

    pub fn merge_single(&mut self, result: &SingleResult) {
        self.merge_stats(FileSlot::SingleOrFirst, &result.stats);
        if result.extracted {
            self.nb_records_extracted += 1;
        }
    }

    pub fn merge_paired(&mut self, result: &PairedResult) {
        self.merge_stats(FileSlot::SingleOrFirst, &result.stats_1);
        self.merge_stats(FileSlot::Second, &result.stats_2);
        if result.extracted {
            self.nb_records_extracted += 2;
        }
    }

    fn merge_stats(&mut self, file_slot: FileSlot, stats: &RecordStats) {
        self.nb_records_tot += 1;
        self.nb_bases += stats.num_bases;
        match file_slot {
            FileSlot::SingleOrFirst => {
                self.nb_hits_tot.0 += stats.hit_count;
                if stats.distinct_record_hit {
                    self.nb_records_hit.0 += 1;
                }
            }
            FileSlot::Second => {
                self.nb_hits_tot.1 += stats.hit_count;
                if stats.distinct_record_hit {
                    self.nb_records_hit.1 += 1;
                }
            }
        }
        for (total, count) in self
            .pattern_hit_counts
            .iter_mut()
            .zip(stats.pattern_hit_counts.iter())
        {
            *total += count;
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::pattern_matching::PatternMatcher;
    use std::time::Duration;

    fn bndmq_processor(logging_active: bool, suppress_output: bool) -> ExtractProcessor {
        let patterns = vec!["ACG".to_string()];
        let matcher = PatternMatcher::new(&patterns, false, false, None).unwrap();
        ExtractProcessor::new(
            Arc::new(matcher),
            patterns.len(),
            logging_active,
            suppress_output,
            false,
        )
    }

    #[test]
    fn single_non_logging_result_uses_first_hit_path_without_hits() {
        let processor = bndmq_processor(false, false);
        let result = processor.process_single_record(
            7,
            RecordInput {
                id: b"seq1",
                seq: b"TTACGTT",
                qual: None,
                format: FastxFormat::Fasta,
            },
        );

        assert_eq!(result.record_index, 7);
        assert!(result.matched);
        assert!(result.extracted);
        assert!(result.hits.is_empty());
        assert_eq!(result.stats.hit_count, 0);
        assert!(!result.stats.distinct_record_hit);
        assert_eq!(result.stats.num_bases, 7);
        assert_eq!(result.output.unwrap().seq, b"TTACGTT");
    }

    #[test]
    fn single_logging_result_reports_all_hits_and_counts() {
        let processor = bndmq_processor(true, false);
        let result = processor.process_single_record(
            0,
            RecordInput {
                id: b"seq1",
                seq: b"ACGACG",
                qual: None,
                format: FastxFormat::Fasta,
            },
        );

        assert!(result.matched);
        assert!(result.extracted);
        assert_eq!(result.stats.hit_count, 2);
        assert!(result.stats.distinct_record_hit);
        assert_eq!(result.stats.pattern_hit_counts, vec![2]);
        assert_eq!(
            result.hits,
            vec![
                RecordHit {
                    file_slot: FileSlot::SingleOrFirst,
                    record_id: b"seq1".to_vec(),
                    pattern_index: 0,
                    position: 0,
                },
                RecordHit {
                    file_slot: FileSlot::SingleOrFirst,
                    record_id: b"seq1".to_vec(),
                    pattern_index: 0,
                    position: 3,
                },
            ]
        );
    }

    #[test]
    fn paired_non_logging_short_circuits_second_read_matching() {
        let processor = bndmq_processor(false, true);
        let result = processor.process_paired_record(
            3,
            RecordInput {
                id: b"seq1/1",
                seq: b"ACG",
                qual: Some(b"III"),
                format: FastxFormat::Fastq,
            },
            RecordInput {
                id: b"seq1/2",
                seq: b"ACG",
                qual: Some(b"III"),
                format: FastxFormat::Fastq,
            },
        );

        assert_eq!(result.pair_index, 3);
        assert!(result.matched);
        assert!(result.extracted);
        assert_eq!(result.stats_1.num_bases, 3);
        assert_eq!(result.stats_2.num_bases, 3);
        assert!(result.hits.is_empty());
        assert!(result.output_1.is_none());
        assert!(result.output_2.is_none());
    }

    #[test]
    fn ordered_result_buffer_drains_contiguous_results_only() {
        let processor = bndmq_processor(false, true);
        let mut buffer = OrderedResultBuffer::new();

        for index in [2, 0] {
            buffer.push(processor.process_single_record(
                index,
                RecordInput {
                    id: b"seq",
                    seq: b"ACG",
                    qual: None,
                    format: FastxFormat::Fasta,
                },
            ));
        }

        let ready = buffer.drain_ready();
        assert_eq!(
            ready
                .iter()
                .map(|result| result.record_index)
                .collect::<Vec<_>>(),
            vec![0]
        );
        assert_eq!(buffer.next_expected_index(), 1);
        assert_eq!(buffer.pending_len(), 1);

        buffer.push(processor.process_single_record(
            1,
            RecordInput {
                id: b"seq",
                seq: b"ACG",
                qual: None,
                format: FastxFormat::Fasta,
            },
        ));

        let ready = buffer.drain_ready();
        assert_eq!(
            ready
                .iter()
                .map(|result| result.record_index)
                .collect::<Vec<_>>(),
            vec![1, 2]
        );
        assert_eq!(buffer.next_expected_index(), 3);
        assert_eq!(buffer.pending_len(), 0);
    }

    #[test]
    fn extract_summary_merges_ordered_single_results() {
        let processor = bndmq_processor(true, false);
        let mut summary = ExtractSummary::new(1);
        let result = processor.process_single_record(
            0,
            RecordInput {
                id: b"seq1",
                seq: b"ACGACG",
                qual: None,
                format: FastxFormat::Fasta,
            },
        );

        summary.merge_single(&result);

        assert_eq!(summary.nb_records_tot, 1);
        assert_eq!(summary.nb_bases, 6);
        assert_eq!(summary.nb_hits_tot, (2, 0));
        assert_eq!(summary.nb_records_hit, (1, 0));
        assert_eq!(summary.nb_records_extracted, 1);
        assert_eq!(summary.pattern_hit_counts, vec![2]);
    }

    #[test]
    fn extract_summary_merges_paired_results_by_file_slot() {
        let processor = bndmq_processor(true, false);
        let mut summary = ExtractSummary::new(1);
        let result = processor.process_paired_record(
            0,
            RecordInput {
                id: b"seq1/1",
                seq: b"ACG",
                qual: Some(b"III"),
                format: FastxFormat::Fastq,
            },
            RecordInput {
                id: b"seq1/2",
                seq: b"TTT",
                qual: Some(b"III"),
                format: FastxFormat::Fastq,
            },
        );

        summary.merge_paired(&result);

        assert_eq!(summary.nb_records_tot, 2);
        assert_eq!(summary.nb_bases, 6);
        assert_eq!(summary.nb_hits_tot, (1, 0));
        assert_eq!(summary.nb_records_hit, (1, 0));
        assert_eq!(summary.nb_records_extracted, 2);
        assert_eq!(summary.pattern_hit_counts, vec![1]);
    }

    #[test]
    fn bounded_pipeline_consumes_out_of_order_results_in_order() {
        let processor = bndmq_processor(false, true);
        let work_items = vec![0_u64, 1, 2, 3];
        let config = PipelineConfig {
            worker_count: 3,
            work_queue_bound: 2,
            result_queue_bound: 2,
        };
        let mut consumed = Vec::new();

        run_bounded_ordered_pipeline(
            work_items,
            config,
            move |index| {
                if index == 0 {
                    thread::sleep(Duration::from_millis(30));
                }
                Ok(vec![processor.process_single_record(
                    index,
                    RecordInput {
                        id: b"seq",
                        seq: b"ACG",
                        qual: None,
                        format: FastxFormat::Fasta,
                    },
                )])
            },
            |result| {
                consumed.push(result.record_index);
                Ok(())
            },
        )
        .unwrap();

        assert_eq!(consumed, vec![0, 1, 2, 3]);
    }

    #[test]
    fn bounded_pipeline_rejects_zero_workers() {
        let error = run_bounded_ordered_pipeline::<u64, SingleResult, _, _>(
            Vec::new(),
            PipelineConfig {
                worker_count: 0,
                work_queue_bound: 1,
                result_queue_bound: 1,
            },
            |_| unreachable!(),
            |_| Ok(()),
        )
        .unwrap_err();

        assert_eq!(
            error.to_string(),
            "Pipeline worker count must be greater than zero."
        );
    }
}
