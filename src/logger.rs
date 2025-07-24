//! Logger utilities

use core::fmt::Write as FmtWrite;
use serde_json::json;
use std::{
    io::{self, Write},
    str,
};

/// A buffered logger that accumulates log records and writes them in batches.
pub struct BufferedLogger {
    buffer: String,
    writer: Option<Box<dyn io::Write>>,
    buffer_size: usize,
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
        self.records.push(record.to_string());
        self.buffer.push_str(record);

        if self.buffer.len() >= self.buffer_size {
            self.flush();
        }
    }

    /// Logs the given fields directly to the buffer without constructing an
    /// intermediate `String` for output.
    pub fn log_fields(&mut self, prefix: &str, record: &[u8], pattern: &str, index: usize) {
        let id_str = str::from_utf8(record).expect("Error during id parsing.");

        // Store the record string for later retrieval
        self.records
            .push(format!("{prefix}\t{id_str}\t{pattern}\t{index}\n"));

        self.buffer.push_str(prefix);
        self.buffer.push('\t');
        self.buffer.push_str(id_str);
        self.buffer.push('\t');
        self.buffer.push_str(pattern);
        self.buffer.push('\t');
        write!(self.buffer, "{index}").unwrap();
        self.buffer.push('\n');

        if self.buffer.len() >= self.buffer_size {
            self.flush();
        }
    }

    /// Writes a header directly to the output without buffering.
    pub fn write_header(&mut self, header: &str) {
        if let Some(writer) = &mut self.writer {
            let _ = writer.write_all(header.as_bytes());
        }
    }

    /// Flushes the buffer to the output.
    pub fn flush(&mut self) {
        if let Some(writer) = &mut self.writer
            && !self.buffer.is_empty()
        {
            let _ = writer.write_all(self.buffer.as_bytes());
            self.buffer.clear();
        }
    }

    /// Returns a reference to the collected records.
    pub fn records(&self) -> &[String] {
        &self.records
    }
}

/// A logger that streams matching records directly to a JSON file.
pub struct JsonLogger {
    buffer: String,
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
            buffer: String::with_capacity(buffer_size),
            writer,
            buffer_size,
            first: true,
        }
    }

    /// Log record fields as a JSON object.
    pub fn log_fields(&mut self, file: &str, record: &[u8], pattern: &str, index: usize) {
        let id_str = str::from_utf8(record).expect("Error during id parsing.");

        if !self.first {
            self.buffer.push_str(",\n");
        }
        self.first = false;

        let value = json!({
            "file": file,
            "record_id": id_str,
            "pattern": pattern,
            "position": index.to_string(),
        });

        let pretty = serde_json::to_string_pretty(&value).unwrap();
        for line in pretty.lines() {
            self.buffer.push_str("    ");
            self.buffer.push_str(line);
            self.buffer.push('\n');
        }

        if self.buffer.len() >= self.buffer_size {
            self.flush();
        }
    }

    /// Flush the internal buffer.
    pub fn flush(&mut self) {
        if let Some(writer) = &mut self.writer
            && !self.buffer.is_empty()
        {
            let _ = writer.write_all(self.buffer.as_bytes());
            self.buffer.clear();
        }
    }

    fn write_indented_value(&mut self, value: &serde_json::Value, indent: usize) {
        let indent_str = " ".repeat(indent);
        let pretty = serde_json::to_string_pretty(value).unwrap();
        for (i, line) in pretty.lines().enumerate() {
            if i > 0 {
                self.buffer.push_str(&indent_str);
            }
            self.buffer.push_str(line);
            self.buffer.push('\n');
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
        self.buffer.push_str("  ],\n  \"meta_information\": ");
        self.write_indented_value(meta_information, 2);
        if self.buffer.ends_with('\n') {
            self.buffer.pop();
        }
        if let Some(stats) = paired_end_stats {
            self.buffer
                .push_str(",\n  \"paired_end_reads_statistics\": ");
            self.write_indented_value(stats, 2);
            if self.buffer.ends_with('\n') {
                self.buffer.pop();
            }
        }
        self.buffer.push_str(",\n  \"pattern_hit_counts\": ");
        self.write_indented_value(pattern_hit_counts, 2);
        if self.buffer.ends_with('\n') {
            self.buffer.pop();
        }
        self.buffer.push_str(",\n  \"summary_statistics\": ");
        self.write_indented_value(summary_statistics, 2);
        if self.buffer.ends_with('\n') {
            self.buffer.pop();
        }
        self.buffer.push_str("\n}\n");
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
}
