#!/usr/bin/env python3

import unittest

import analyze


class LowerEnvelopeTests(unittest.TestCase):
    def test_single_algorithm_covers_all_corpus_sizes(self):
        regions = analyze.lower_envelope({"hash": (10.0, 2.0)})

        self.assertEqual(
            regions,
            [{"algorithm": "hash", "start": 0.0, "end": None}],
        )

    def test_fixed_cost_winner_yields_to_faster_search(self):
        regions = analyze.lower_envelope(
            {
                "hash": (10.0, 2.0),
                "aho_corasick": (110.0, 1.0),
            }
        )

        self.assertEqual(
            [region["algorithm"] for region in regions],
            ["hash", "aho_corasick"],
        )
        self.assertAlmostEqual(regions[0]["end"], 100.0)
        self.assertAlmostEqual(regions[1]["start"], 100.0)
        self.assertIsNone(regions[1]["end"])

    def test_dominated_line_is_excluded(self):
        regions = analyze.lower_envelope(
            {
                "hash": (10.0, 2.0),
                "bndmq": (20.0, 2.5),
                "aho_corasick": (110.0, 1.0),
            }
        )

        self.assertEqual(
            [region["algorithm"] for region in regions],
            ["hash", "aho_corasick"],
        )

    def test_three_algorithms_can_form_two_transitions(self):
        regions = analyze.lower_envelope(
            {
                "hash": (0.0, 4.0),
                "bndmq": (100.0, 2.0),
                "aho_corasick": (400.0, 1.0),
            }
        )

        self.assertEqual(
            [region["algorithm"] for region in regions],
            ["hash", "bndmq", "aho_corasick"],
        )
        self.assertAlmostEqual(regions[0]["end"], 50.0)
        self.assertAlmostEqual(regions[1]["end"], 300.0)


if __name__ == "__main__":
    unittest.main()
