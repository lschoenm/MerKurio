//! Logger utilities

use core::fmt::Write as FmtWrite;
use std::{
    io::{self, Write},
    str,
};

pub fn append_log_fields(
    buffer: &mut String,
    prefix: &str,
    record: &[u8],
    pattern: &str,
    index: usize,
) {
    let id_str = str::from_utf8(record).expect("Error during id parsing.");
    buffer.push_str(prefix);
    buffer.push('\t');
    buffer.push_str(id_str);
    buffer.push('\t');
    buffer.push_str(pattern);
    buffer.push('\t');
    write!(buffer, "{index}").unwrap();
    buffer.push('\n');
}

pub fn append_json_log_fields(
    buffer: &mut Vec<u8>,
    first: &mut bool,
    file: &str,
    record: &[u8],
    pattern: &str,
    index: usize,
) {
    let id_str = str::from_utf8(record).expect("Error during id parsing.");

    if !*first {
        buffer.extend_from_slice(b",\n");
    }
    *first = false;

    buffer.extend_from_slice(b"    {\n      \"file\": ");
    serde_json::to_writer(&mut *buffer, file).unwrap();
    buffer.extend_from_slice(b",\n      \"pattern\": ");
    serde_json::to_writer(&mut *buffer, pattern).unwrap();
    buffer.extend_from_slice(b",\n      \"position\": \"");
    write!(buffer, "{index}").unwrap();
    buffer.extend_from_slice(b"\",\n      \"record_id\": ");
    serde_json::to_writer(&mut *buffer, id_str).unwrap();
    buffer.extend_from_slice(b"\n    }\n");
}

/// A buffered logger that streams log records in batches.
pub struct BufferedLogger {
    buffer: String,
    writer: Option<Box<dyn io::Write>>,
    buffer_size: usize,
    /// Retained only for records-only/test mode when no writer is configured.
    records: Vec<String>,
}

impl BufferedLogger {
    /// Creates a new `BufferedLogger` with the given writer and buffer size.
    pub fn new(writer: Option<Box<dyn io::Write>>, buffer_size: usize) -> Self {
        Self {
            buffer: String::with_capacity(buffer_size),
            writer,
            buffer_size,
            records: Vec::new(),
        }
    }

    /// Logs a record to the buffer and writes to output if buffer is full.
    pub fn log_record(&mut self, record: &str) {
        if self.writer.is_none() {
            self.records.push(record.to_string());
            return;
        }

        self.buffer.push_str(record);

        if self.buffer.len() >= self.buffer_size {
            self.flush();
        }
    }

    /// Logs the given fields directly to the buffer without constructing an
    /// intermediate `String` for output.
    pub fn log_fields(&mut self, prefix: &str, record: &[u8], pattern: &str, index: usize) {
        if self.writer.is_none() {
            let id_str = str::from_utf8(record).expect("Error during id parsing.");
            self.records
                .push(format!("{prefix}\t{id_str}\t{pattern}\t{index}\n"));
            return;
        }

        append_log_fields(&mut self.buffer, prefix, record, pattern, index);

        if self.buffer.len() >= self.buffer_size {
            self.flush();
        }
    }

    pub fn log_fragment(&mut self, fragment: &str) -> io::Result<()> {
        if fragment.is_empty() {
            return Ok(());
        }
        if self.writer.is_none() {
            self.records.push(fragment.to_string());
            return Ok(());
        }
        if fragment.len() >= self.buffer_size {
            self.try_flush()?;
            if let Some(writer) = &mut self.writer {
                writer.write_all(fragment.as_bytes())?;
            }
            return Ok(());
        }
        self.buffer.push_str(fragment);
        if self.buffer.len() >= self.buffer_size {
            self.try_flush()?;
        }
        Ok(())
    }

    /// Writes a header directly to the output without buffering.
    pub fn write_header(&mut self, header: &str) {
        if let Some(writer) = &mut self.writer {
            let _ = writer.write_all(header.as_bytes());
        }
    }

    /// Flushes the buffer to the output.
    pub fn flush(&mut self) {
        let _ = self.try_flush();
    }

    fn try_flush(&mut self) -> io::Result<()> {
        if let Some(writer) = &mut self.writer
            && !self.buffer.is_empty()
        {
            writer.write_all(self.buffer.as_bytes())?;
            self.buffer.clear();
        }
        Ok(())
    }

    /// Returns records collected in records-only mode.
    pub fn records(&self) -> &[String] {
        &self.records
    }
}

/// A logger that streams matching records directly to a JSON file.
pub struct JsonLogger {
    buffer: Vec<u8>,
    writer: Option<Box<dyn io::Write>>,
    buffer_size: usize,
    first: bool,
}

impl JsonLogger {
    /// Create a new `JsonLogger` with the given writer and buffer size.
    pub fn new(mut writer: Option<Box<dyn io::Write>>, buffer_size: usize) -> Self {
        if let Some(w) = &mut writer {
            let _ = w.write_all(b"{\n  \"matching_records\": [\n");
        }
        Self {
            buffer: Vec::with_capacity(buffer_size),
            writer,
            buffer_size,
            first: true,
        }
    }

    /// Log record fields as a JSON object.
    pub fn log_fields(&mut self, file: &str, record: &[u8], pattern: &str, index: usize) {
        append_json_log_fields(
            &mut self.buffer,
            &mut self.first,
            file,
            record,
            pattern,
            index,
        );

        if self.buffer.len() >= self.buffer_size {
            self.flush();
        }
    }

    pub fn log_fragment(&mut self, fragment: &[u8]) -> io::Result<()> {
        if fragment.is_empty() {
            return Ok(());
        }
        if !self.first {
            self.buffer.extend_from_slice(b",\n");
        }
        self.first = false;
        if fragment.len() >= self.buffer_size {
            self.try_flush()?;
            if let Some(writer) = &mut self.writer {
                writer.write_all(fragment)?;
            }
            return Ok(());
        }
        self.buffer.extend_from_slice(fragment);
        if self.buffer.len() >= self.buffer_size {
            self.try_flush()?;
        }
        Ok(())
    }

    /// Flush the internal buffer.
    pub fn flush(&mut self) {
        let _ = self.try_flush();
    }

    fn try_flush(&mut self) -> io::Result<()> {
        if let Some(writer) = &mut self.writer
            && !self.buffer.is_empty()
        {
            writer.write_all(&self.buffer)?;
            self.buffer.clear();
        }
        Ok(())
    }

    fn write_indented_value(&mut self, value: &serde_json::Value, indent: usize) {
        let indent_str = " ".repeat(indent);
        let pretty = serde_json::to_string_pretty(value).unwrap();
        for (i, line) in pretty.lines().enumerate() {
            if i > 0 {
                self.buffer.extend_from_slice(indent_str.as_bytes());
            }
            self.buffer.extend_from_slice(line.as_bytes());
            self.buffer.push(b'\n');
        }
    }

    /// Finalize the JSON output by writing summary information.
    pub fn finalize(
        mut self,
        meta_information: &serde_json::Value,
        pattern_hit_counts: &serde_json::Value,
        summary_statistics: &serde_json::Value,
        paired_end_stats: Option<&serde_json::Value>,
    ) {
        self.buffer
            .extend_from_slice(b"  ],\n  \"meta_information\": ");
        self.write_indented_value(meta_information, 2);
        if self.buffer.ends_with(b"\n") {
            self.buffer.pop();
        }
        if let Some(stats) = paired_end_stats {
            self.buffer
                .extend_from_slice(b",\n  \"paired_end_reads_statistics\": ");
            self.write_indented_value(stats, 2);
            if self.buffer.ends_with(b"\n") {
                self.buffer.pop();
            }
        }
        self.buffer
            .extend_from_slice(b",\n  \"pattern_hit_counts\": ");
        self.write_indented_value(pattern_hit_counts, 2);
        if self.buffer.ends_with(b"\n") {
            self.buffer.pop();
        }
        self.buffer
            .extend_from_slice(b",\n  \"summary_statistics\": ");
        self.write_indented_value(summary_statistics, 2);
        if self.buffer.ends_with(b"\n") {
            self.buffer.pop();
        }
        self.buffer.extend_from_slice(b"\n}\n");
        self.flush();
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_buffered_logger_basic() {
        let mut logger = BufferedLogger::new(None, 1024);

        // Test header writing
        logger.write_header("Header 1\n");
        logger.write_header("Header 2\n");
        logger.flush();

        // Test record logging
        logger.log_record("Record 1\n");
        logger.log_record("Record 2\n");
        logger.flush();

        // Verify records are stored correctly
        let records = logger.records();
        assert_eq!(records.len(), 2);
        assert_eq!(records[0], "Record 1\n");
        assert_eq!(records[1], "Record 2\n");
    }

    #[test]
    fn test_buffered_logger_no_writer() {
        let mut logger = BufferedLogger::new(None, 1024);

        // These should not panic
        logger.write_header("Header\n");
        logger.log_record("Record\n");
        logger.flush();

        // Records should still be stored
        let records = logger.records();
        assert_eq!(records.len(), 1);
        assert_eq!(records[0], "Record\n");
    }

    #[test]
    fn test_buffered_logger_empty() {
        let mut logger = BufferedLogger::new(None, 1024);

        // Flush empty logger
        logger.flush();
        assert_eq!(logger.records().len(), 0);
    }

    #[test]
    fn test_buffered_logger_records_only() {
        let mut logger = BufferedLogger::new(None, 1024);

        // Log some records
        logger.log_record("Record 1\n");
        logger.log_record("Record 2\n");
        logger.log_record("Record 3\n");

        // Check if records are stored correctly
        let records = logger.records();
        assert_eq!(records.len(), 3);
        assert_eq!(records[0], "Record 1\n");
        assert_eq!(records[1], "Record 2\n");
        assert_eq!(records[2], "Record 3\n");
    }

    #[test]
    fn test_buffered_logger_with_writer_does_not_retain_records() {
        let mut logger = BufferedLogger::new(Some(Box::new(std::io::sink())), 1024);

        logger.log_record("Record 1\n");
        logger.log_fields("file.fastq", b"read-1", "ACGT", 7);
        logger.flush();

        assert!(logger.records().is_empty());
    }
}
