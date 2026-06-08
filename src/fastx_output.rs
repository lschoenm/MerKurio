use anyhow::Result;
use std::io::Write;

#[derive(Debug, Clone, Copy, Eq, PartialEq)]
pub enum FastxFormat {
    Fasta,
    Fastq,
}

#[derive(Debug, Clone, Copy)]
pub struct FastxRecordView<'a> {
    pub id: &'a [u8],
    pub seq: &'a [u8],
    pub qual: Option<&'a [u8]>,
    pub format: FastxFormat,
}

/// Write a normalized FASTA/FASTQ record.
///
/// This intentionally does not preserve original line wrapping or line endings.
/// It is the output layer for parser-independent extract implementations.
pub fn write_fastx_record(writer: &mut dyn Write, record: FastxRecordView<'_>) -> Result<()> {
    match record.format {
        FastxFormat::Fasta => write_fasta_record(writer, record.id, record.seq),
        FastxFormat::Fastq => {
            let qual = record
                .qual
                .ok_or_else(|| anyhow::anyhow!("FASTQ output requires quality scores."))?;
            if qual.len() != record.seq.len() {
                anyhow::bail!(
                    "FASTQ sequence and quality lengths differ: {} sequence bases, {} quality scores.",
                    record.seq.len(),
                    qual.len()
                );
            }
            write_fastq_record(writer, record.id, record.seq, qual)
        }
    }
}

fn write_fasta_record(writer: &mut dyn Write, id: &[u8], seq: &[u8]) -> Result<()> {
    writer.write_all(b">")?;
    writer.write_all(id)?;
    writer.write_all(b"\n")?;
    writer.write_all(seq)?;
    writer.write_all(b"\n")?;
    Ok(())
}

fn write_fastq_record(writer: &mut dyn Write, id: &[u8], seq: &[u8], qual: &[u8]) -> Result<()> {
    writer.write_all(b"@")?;
    writer.write_all(id)?;
    writer.write_all(b"\n")?;
    writer.write_all(seq)?;
    writer.write_all(b"\n+\n")?;
    writer.write_all(qual)?;
    writer.write_all(b"\n")?;
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn write_normalized_fasta_record() {
        let mut output = Vec::new();

        write_fastx_record(
            &mut output,
            FastxRecordView {
                id: b"seq1 description",
                seq: b"ACGTACGT",
                qual: None,
                format: FastxFormat::Fasta,
            },
        )
        .unwrap();

        assert_eq!(output, b">seq1 description\nACGTACGT\n");
    }

    #[test]
    fn write_normalized_fastq_record() {
        let mut output = Vec::new();

        write_fastx_record(
            &mut output,
            FastxRecordView {
                id: b"seq1/1",
                seq: b"ACTTACGT",
                qual: Some(b"IIIIIIII"),
                format: FastxFormat::Fastq,
            },
        )
        .unwrap();

        assert_eq!(output, b"@seq1/1\nACTTACGT\n+\nIIIIIIII\n");
    }

    #[test]
    fn fastq_requires_quality_scores() {
        let mut output = Vec::new();

        let error = write_fastx_record(
            &mut output,
            FastxRecordView {
                id: b"seq1/1",
                seq: b"ACTTACGT",
                qual: None,
                format: FastxFormat::Fastq,
            },
        )
        .unwrap_err();

        assert_eq!(error.to_string(), "FASTQ output requires quality scores.");
        assert!(output.is_empty());
    }

    #[test]
    fn fastq_rejects_quality_length_mismatch() {
        let mut output = Vec::new();

        let error = write_fastx_record(
            &mut output,
            FastxRecordView {
                id: b"seq1/1",
                seq: b"ACTTACGT",
                qual: Some(b"IIII"),
                format: FastxFormat::Fastq,
            },
        )
        .unwrap_err();

        assert_eq!(
            error.to_string(),
            "FASTQ sequence and quality lengths differ: 8 sequence bases, 4 quality scores."
        );
        assert!(output.is_empty());
    }
}
