#!/usr/bin/env python3

import unittest

import refine


class RefinementTargetTests(unittest.TestCase):
    def test_geometric_midpoint_stays_inside_integer_bracket(self):
        self.assertEqual(refine.geometric_midpoint(16, 32), 23)
        self.assertEqual(refine.geometric_midpoint(2, 4), 3)

    def test_transition_targets_only_changed_winners(self):
        grouped = {
            (31, "first"): [
                {"patterns": 8, "algorithm": "bndmq"},
                {"patterns": 16, "algorithm": "bndmq"},
                {"patterns": 32, "algorithm": "aho_corasick"},
                {"patterns": 64, "algorithm": "aho_corasick"},
            ]
        }

        targets = refine.transition_targets(grouped)

        self.assertEqual(
            targets,
            [
                {
                    "axis": "patterns",
                    "k": 31,
                    "mode": "first",
                    "patterns": 23,
                    "lower_k": 31,
                    "upper_k": 31,
                    "lower_patterns": 16,
                    "upper_patterns": 32,
                    "lower_algorithm": "bndmq",
                    "upper_algorithm": "aho_corasick",
                }
            ],
        )

    def test_adjacent_integer_transition_cannot_be_refined(self):
        grouped = {
            (80, "all"): [
                {"patterns": 1, "algorithm": "aho_corasick"},
                {"patterns": 2, "algorithm": "hash"},
            ]
        }

        self.assertEqual(refine.transition_targets(grouped), [])

    def test_pattern_length_transition_uses_arithmetic_midpoint(self):
        grouped = {
            (24, "all"): [{"patterns": 16, "algorithm": "aho_corasick"}],
            (31, "all"): [{"patterns": 16, "algorithm": "bndmq"}],
        }

        targets = refine.transition_targets(grouped)

        self.assertEqual(len(targets), 1)
        self.assertEqual(targets[0]["axis"], "k")
        self.assertEqual(targets[0]["k"], 27)
        self.assertEqual(targets[0]["patterns"], 16)

    def test_pattern_length_transition_anchors_bndmq_boundary_at_65(self):
        grouped = {
            (64, "first"): [{"patterns": 32, "algorithm": "bndmq"}],
            (80, "first"): [{"patterns": 32, "algorithm": "hash"}],
        }

        targets = refine.transition_targets(grouped)

        self.assertEqual(targets[0]["k"], 65)


if __name__ == "__main__":
    unittest.main()
