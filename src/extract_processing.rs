use crossbeam_channel::{Sender, bounded};
use std::collections::BTreeMap;
use std::sync::Arc;
use std::thread;

#[derive(Debug, Clone, Copy, Eq, PartialEq)]
pub enum FileSlot {
    SingleOrFirst,
    Second,
}

pub trait IndexedResult {
    fn index(&self) -> u64;

    fn index_span(&self) -> u64 {
        1
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

    pub fn pop_ready(&mut self) -> Option<T> {
        let result = self.pending.remove(&self.next_expected_index)?;
        self.next_expected_index += result.index_span();
        Some(result)
    }
}

impl<T> Default for OrderedResultBuffer<T> {
    fn default() -> Self {
        Self::new()
    }
}

#[derive(Debug)]
struct ResultChunk<R> {
    start_index: u64,
    results: Vec<R>,
}

impl<R: IndexedResult> IndexedResult for ResultChunk<R> {
    fn index(&self) -> u64 {
        self.start_index
    }

    fn index_span(&self) -> u64 {
        self.results.iter().map(IndexedResult::index_span).sum()
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
    let (result_tx, result_rx) = bounded::<ResultChunk<R>>(config.result_queue_bound);

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
                let results = process_work(work)?;
                if let Some(first_result) = results.first() {
                    result_tx
                        .send(ResultChunk {
                            start_index: first_result.index(),
                            results,
                        })
                        .map_err(|_| {
                            anyhow::anyhow!(
                                "Pipeline result queue closed before all result chunks were sent."
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
    while let Ok(result_chunk) = result_rx.recv() {
        ordered_results.push(result_chunk);
        while let Some(ready_chunk) = ordered_results.pop_ready() {
            for ready in ready_chunk.results {
                if let Err(error) = consume_ready(ready) {
                    consume_error = Some(error);
                    break;
                }
            }
            if consume_error.is_some() {
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
            "Pipeline completed with {} out-of-order result chunk(s) still pending; missing result index {}.",
            ordered_results.pending_len(),
            ordered_results.next_expected_index()
        );
    }

    Ok(())
}

pub fn run_bounded_ordered_pipeline_with_producer<W, R, PR, P, C>(
    config: PipelineConfig,
    produce_work: PR,
    process_work: P,
    mut consume_ready: C,
) -> anyhow::Result<()>
where
    W: Send + 'static,
    R: IndexedResult + Send + 'static,
    PR: FnOnce(Sender<W>) -> anyhow::Result<()> + Send + 'static,
    P: Fn(W) -> anyhow::Result<R> + Send + Sync + 'static,
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

    let producer = thread::spawn(move || produce_work(work_tx));

    let process_work = Arc::new(process_work);
    let mut workers = Vec::with_capacity(config.worker_count);
    for _ in 0..config.worker_count {
        let work_rx = work_rx.clone();
        let result_tx = result_tx.clone();
        let process_work = Arc::clone(&process_work);
        workers.push(thread::spawn(move || -> anyhow::Result<()> {
            while let Ok(work) = work_rx.recv() {
                let result = process_work(work)?;
                result_tx.send(result).map_err(|_| {
                    anyhow::anyhow!(
                        "Pipeline result queue closed before all result chunks were sent."
                    )
                })?;
            }
            Ok(())
        }));
    }
    drop(result_tx);

    let mut ordered_results = OrderedResultBuffer::new();
    let mut consume_error = None;
    while let Ok(result) = result_rx.recv() {
        ordered_results.push(result);
        while let Some(ready) = ordered_results.pop_ready() {
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
            "Pipeline completed with {} out-of-order result chunk(s) still pending; missing result index {}.",
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

    pub fn record_searched(&mut self, num_bases: usize) {
        self.nb_records_tot += 1;
        self.nb_bases += num_bases;
    }

    pub fn record_hit(&mut self, file_slot: FileSlot) {
        match file_slot {
            FileSlot::SingleOrFirst => self.nb_records_hit.0 += 1,
            FileSlot::Second => self.nb_records_hit.1 += 1,
        }
    }

    pub fn pattern_hit(&mut self, file_slot: FileSlot, pattern_index: usize) {
        match file_slot {
            FileSlot::SingleOrFirst => self.nb_hits_tot.0 += 1,
            FileSlot::Second => self.nb_hits_tot.1 += 1,
        }
        self.pattern_hit_counts[pattern_index] += 1;
    }

    pub fn extracted_records(&mut self, count: usize) {
        self.nb_records_extracted += count;
    }

    pub fn merge(&mut self, other: &Self) {
        self.nb_records_tot += other.nb_records_tot;
        self.nb_bases += other.nb_bases;
        self.nb_hits_tot.0 += other.nb_hits_tot.0;
        self.nb_hits_tot.1 += other.nb_hits_tot.1;
        self.nb_records_hit.0 += other.nb_records_hit.0;
        self.nb_records_hit.1 += other.nb_records_hit.1;
        self.nb_records_extracted += other.nb_records_extracted;
        debug_assert_eq!(
            self.pattern_hit_counts.len(),
            other.pattern_hit_counts.len()
        );
        for (total, chunk_count) in self
            .pattern_hit_counts
            .iter_mut()
            .zip(&other.pattern_hit_counts)
        {
            *total += chunk_count;
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::time::Duration;

    #[derive(Debug)]
    struct TestResult {
        index: u64,
        span: u64,
    }

    impl IndexedResult for TestResult {
        fn index(&self) -> u64 {
            self.index
        }

        fn index_span(&self) -> u64 {
            self.span
        }
    }

    #[test]
    fn ordered_result_buffer_pops_contiguous_results_only() {
        let mut buffer = OrderedResultBuffer::new();

        for index in [2, 0] {
            buffer.push(TestResult { index, span: 1 });
        }

        assert_eq!(buffer.pop_ready().map(|result| result.index), Some(0));
        assert!(buffer.pop_ready().is_none());
        assert_eq!(buffer.next_expected_index(), 1);
        assert_eq!(buffer.pending_len(), 1);

        buffer.push(TestResult { index: 1, span: 1 });

        assert_eq!(buffer.pop_ready().map(|result| result.index), Some(1));
        assert_eq!(buffer.pop_ready().map(|result| result.index), Some(2));
        assert!(buffer.pop_ready().is_none());
        assert_eq!(buffer.next_expected_index(), 3);
        assert_eq!(buffer.pending_len(), 0);
    }

    #[test]
    fn extract_summary_merges_chunk_totals() {
        let mut total = ExtractSummary::new(2);
        let mut chunk = ExtractSummary::new(2);
        chunk.record_searched(6);
        chunk.record_searched(4);
        chunk.pattern_hit(FileSlot::SingleOrFirst, 0);
        chunk.pattern_hit(FileSlot::SingleOrFirst, 0);
        chunk.pattern_hit(FileSlot::Second, 1);
        chunk.record_hit(FileSlot::SingleOrFirst);
        chunk.record_hit(FileSlot::Second);
        chunk.extracted_records(2);

        total.merge(&chunk);

        assert_eq!(total.nb_records_tot, 2);
        assert_eq!(total.nb_bases, 10);
        assert_eq!(total.nb_hits_tot, (2, 1));
        assert_eq!(total.nb_records_hit, (1, 1));
        assert_eq!(total.nb_records_extracted, 2);
        assert_eq!(total.pattern_hit_counts, vec![2, 1]);
    }

    #[test]
    fn bounded_pipeline_consumes_out_of_order_results_in_order() {
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
                Ok(vec![TestResult { index, span: 1 }])
            },
            |result| {
                consumed.push(result.index);
                Ok(())
            },
        )
        .unwrap();

        assert_eq!(consumed, vec![0, 1, 2, 3]);
    }

    #[test]
    fn bounded_pipeline_consumes_chunked_results_in_order() {
        let work_items = vec![2_u64, 0];
        let config = PipelineConfig {
            worker_count: 2,
            work_queue_bound: 2,
            result_queue_bound: 2,
        };
        let mut consumed = Vec::new();

        run_bounded_ordered_pipeline(
            work_items,
            config,
            move |start_index| {
                if start_index == 0 {
                    thread::sleep(Duration::from_millis(30));
                }
                Ok((start_index..start_index + 2)
                    .map(|index| TestResult { index, span: 1 })
                    .collect())
            },
            |result| {
                consumed.push(result.index);
                Ok(())
            },
        )
        .unwrap();

        assert_eq!(consumed, vec![0, 1, 2, 3]);
    }

    #[test]
    fn bounded_pipeline_rejects_zero_workers() {
        let error = run_bounded_ordered_pipeline::<u64, TestResult, _, _>(
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
