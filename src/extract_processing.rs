#[derive(Debug, Clone, Copy, Eq, PartialEq)]
pub enum FileSlot {
    SingleOrFirst,
    Second,
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
}
