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
                    "k": 31,
                    "mode": "first",
                    "patterns": 23,
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


if __name__ == "__main__":
    unittest.main()
