use paraseq::{fastx, Record};
use std::sync::{
    atomic::{AtomicUsize, Ordering},
    Arc,
};

#[derive(Debug, Eq, PartialEq)]
struct ParsedRecord {
    id: Vec<u8>,
    seq: Vec<u8>,
    seq_raw: Vec<u8>,
    qual: Option<Vec<u8>>,
}

fn read_fastx(path: &str) -> Vec<ParsedRecord> {
    let mut reader = fastx::Reader::from_path(path).unwrap();
    let mut record_set = reader.new_record_set();
    let mut records = Vec::new();

    while record_set.fill(&mut reader).unwrap() {
        for record in record_set.iter() {
            let record = record.unwrap();
            records.push(ParsedRecord {
                id: record.id().to_vec(),
                seq: record.seq().into_owned(),
                seq_raw: record.seq_raw().to_vec(),
                qual: record.qual().map(|qual| qual.to_vec()),
            });
        }
    }

    records
}

#[derive(Clone, Default)]
struct CountingPairedProcessor {
    pairs: Arc<AtomicUsize>,
}

impl<R: Record> paraseq::prelude::PairedParallelProcessor<R> for CountingPairedProcessor {
    fn process_record_pair(&mut self, record_1: R, record_2: R) -> paraseq::Result<()> {
        assert!(record_1.id().ends_with(b"/1"));
        assert!(record_2.id().ends_with(b"/2"));
        self.pairs.fetch_add(1, Ordering::Relaxed);
        Ok(())
    }
}

#[test]
fn paraseq_reads_single_end_fasta_fixture() {
    let records = read_fastx("tests/fixtures/input/simple.fasta");

    assert_eq!(records.len(), 3);
    assert_eq!(records[0].id, b"seq1");
    assert_eq!(records[0].seq, b"ACGTACGT");
    assert_eq!(records[0].seq_raw, b"ACGTACGT\n");
    assert_eq!(records[0].qual, None);
}

#[test]
fn paraseq_reads_single_end_fastq_fixture() {
    let records = read_fastx("tests/fixtures/input/paired-1.fastq");

    assert_eq!(records.len(), 2);
    assert_eq!(records[0].id, b"seq1/1");
    assert_eq!(records[0].seq, b"ACTTACGT");
    assert_eq!(records[0].seq_raw, b"ACTTACGT");
    assert_eq!(records[0].qual.as_deref(), Some(&b"IIIIIIII"[..]));
}

#[test]
fn paraseq_strips_multiline_fasta_seq_but_keeps_raw_seq() {
    let records = read_fastx("tests/fixtures/input/fixed-width.faa");

    assert_eq!(records.len(), 1);
    assert_eq!(records[0].id, b"protein1");
    assert_eq!(records[0].seq.len(), 280);
    assert_eq!(records[0].seq_raw.len(), 284);
    assert!(records[0].seq_raw.contains(&b'\n'));
    assert!(!records[0].seq.contains(&b'\n'));
}

#[test]
fn paraseq_reads_existing_compressed_fasta_inputs() {
    let gz = read_fastx("tests/data/sample.fasta.gz");
    let bz2 = read_fastx("tests/data/sample.fasta.bz2");
    let xz = read_fastx("tests/data/sample.fasta.xz");

    assert!(!gz.is_empty());
    assert_eq!(gz, bz2);
    assert_eq!(gz, xz);
}

#[test]
fn paraseq_processes_two_file_paired_end_fastq() {
    let collection = fastx::Collection::from_paths(
        &[
            "tests/fixtures/input/paired-1.fastq",
            "tests/fixtures/input/paired-2.fastq",
        ],
        fastx::CollectionType::Paired,
    )
    .unwrap();
    let mut processor = CountingPairedProcessor::default();

    collection
        .process_parallel_paired(&mut processor, 1, None)
        .unwrap();

    assert_eq!(processor.pairs.load(Ordering::Relaxed), 2);
}

#[test]
fn paraseq_errors_on_mismatched_paired_end_lengths() {
    let temp_dir = tempfile::tempdir().unwrap();
    let read_1 = temp_dir.path().join("read_1.fastq");
    let read_2 = temp_dir.path().join("read_2.fastq");
    std::fs::write(
        &read_1,
        b"@seq1/1\nACTTACGT\n+\nIIIIIIII\n@seq2/1\nTTTTTTTT\n+\nIIIIIIII\n",
    )
    .unwrap();
    std::fs::write(&read_2, b"@seq1/2\nGCTATAAT\n+\nIIIIIIII\n").unwrap();

    let collection = fastx::Collection::from_paths(
        &[read_1.to_string_lossy().as_ref(), read_2.to_string_lossy().as_ref()],
        fastx::CollectionType::Paired,
    )
    .unwrap();
    let mut processor = CountingPairedProcessor::default();

    let error = collection
        .process_parallel_paired(&mut processor, 1, None)
        .unwrap_err();

    assert!(error.to_string().contains("Incompatible record set sizes"));
}
