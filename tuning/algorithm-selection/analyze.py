#!/usr/bin/env python3

import csv
import statistics
from collections import defaultdict
from pathlib import Path


ROOT = Path(__file__).resolve().parent
RESULTS = ROOT / "results"
RAW = RESULTS / "algorithm_sweep.csv"
STATUS = RESULTS / "cell_status.csv"
SUMMARY = RESULTS / "algorithm_summary.csv"
WINNERS = RESULTS / "algorithm_winners.csv"
SELECTION_MAP = RESULTS / "selection_map.svg"
MODES = ("first", "all")
ALGORITHMS = ("bndmq", "hash", "aho_corasick")
COLORS = {
    "bndmq": "#1a73e8",
    "hash": "#188038",
    "aho_corasick": "#d93025",
}
LABELS = {
    "bndmq": "B",
    "hash": "H",
    "aho_corasick": "A",
}


def load_raw():
    timings = defaultdict(list)
    build_times = {}
    validation = defaultdict(dict)
    with RAW.open(newline="") as handle:
        for row in csv.DictReader(handle):
            key = (
                int(row["k"]),
                int(row["patterns"]),
                row["mode"],
                row["algorithm"],
            )
            timings[key].append(float(row["ns_per_base"]))
            build_times[key] = int(row["build_ns"])
            cell = key[:3]
            observed = (
                int(row["validation_checksum"]),
                int(row["matching_records"]),
                int(row["matches"]),
            )
            previous = validation[cell].setdefault(key[3], observed)
            if previous != observed:
                raise SystemExit(f"Inconsistent validation result within {key}")
    return timings, build_times, validation


def load_statuses():
    statuses = {}
    with STATUS.open(newline="") as handle:
        for row in csv.DictReader(handle):
            key = (
                int(row["k"]),
                int(row["patterns"]),
                row["mode"],
                row["algorithm"],
            )
            statuses[key] = row
    return statuses


def validate_algorithms(validation):
    for cell, by_algorithm in validation.items():
        observed = set(by_algorithm.values())
        if len(observed) > 1:
            details = ", ".join(
                f"{algorithm}={value}"
                for algorithm, value in sorted(by_algorithm.items())
            )
            raise SystemExit(f"Algorithm result mismatch at {cell}: {details}")


def summarize(timings, build_times):
    summary = {}
    for key, values in timings.items():
        summary[key] = {
            "runs": len(values),
            "median": statistics.median(values),
            "mean": statistics.fmean(values),
            "stddev": statistics.stdev(values) if len(values) > 1 else 0.0,
            "minimum": min(values),
            "maximum": max(values),
            "build_ns": build_times[key],
        }
    return summary


def write_tables(summary):
    best_by_cell = {}
    for (k, patterns, mode, _algorithm), values in summary.items():
        cell = (k, patterns, mode)
        best_by_cell[cell] = min(
            values["median"], best_by_cell.get(cell, float("inf"))
        )

    with SUMMARY.open("w", newline="") as handle:
        writer = csv.writer(handle)
        writer.writerow(
            [
                "k",
                "patterns",
                "mode",
                "algorithm",
                "runs",
                "build_ms",
                "median_ns_per_base",
                "mean_ns_per_base",
                "stddev_ns_per_base",
                "min_ns_per_base",
                "max_ns_per_base",
                "slower_than_best_pct",
            ]
        )
        for key, values in sorted(summary.items()):
            k, patterns, mode, algorithm = key
            relative = (
                values["median"] / best_by_cell[(k, patterns, mode)] - 1.0
            ) * 100.0
            writer.writerow(
                [
                    k,
                    patterns,
                    mode,
                    algorithm,
                    values["runs"],
                    f"{values['build_ns'] / 1_000_000.0:.3f}",
                    f"{values['median']:.9f}",
                    f"{values['mean']:.9f}",
                    f"{values['stddev']:.9f}",
                    f"{values['minimum']:.9f}",
                    f"{values['maximum']:.9f}",
                    f"{relative:.4f}",
                ]
            )

    winners = {}
    with WINNERS.open("w", newline="") as handle:
        writer = csv.writer(handle)
        writer.writerow(
            [
                "k",
                "patterns",
                "mode",
                "best_algorithm",
                "median_ns_per_base",
                "second_algorithm",
                "margin_pct",
            ]
        )
        for cell in sorted(best_by_cell):
            k, patterns, mode = cell
            choices = sorted(
                (
                    values["median"],
                    algorithm,
                )
                for (
                    candidate_k,
                    candidate_patterns,
                    candidate_mode,
                    algorithm,
                ), values in summary.items()
                if (candidate_k, candidate_patterns, candidate_mode) == cell
            )
            best_value, best_algorithm = choices[0]
            if len(choices) > 1:
                second_value, second_algorithm = choices[1]
                margin = (second_value / best_value - 1.0) * 100.0
            else:
                second_algorithm = ""
                margin = None
            winners[cell] = {
                "algorithm": best_algorithm,
                "margin": margin,
                "successful_algorithms": len(choices),
            }
            writer.writerow(
                [
                    k,
                    patterns,
                    mode,
                    best_algorithm,
                    f"{best_value:.9f}",
                    second_algorithm,
                    "" if margin is None else f"{margin:.4f}",
                ]
            )
    return winners


def write_selection_map(winners, statuses):
    k_values = sorted({key[0] for key in statuses})
    pattern_counts = sorted({key[1] for key in statuses})
    modes = [mode for mode in MODES if any(key[2] == mode for key in statuses)]
    cell_width = max(58, min(92, 1120 // max(1, len(pattern_counts))))
    cell_height = max(25, min(38, 500 // max(1, len(k_values))))
    left = 96
    right = 28
    top = 130
    panel_gap = 100
    panel_height = len(k_values) * cell_height
    width = max(720, left + len(pattern_counts) * cell_width + right)
    height = top + len(modes) * panel_height + max(0, len(modes) - 1) * panel_gap + 100

    svg = [
        f'<svg xmlns="http://www.w3.org/2000/svg" width="{width}" height="{height}" '
        f'viewBox="0 0 {width} {height}">',
        '<rect width="100%" height="100%" fill="#ffffff"/>',
        '<style>text{font-family:-apple-system,BlinkMacSystemFont,"Segoe UI",sans-serif;'
        'fill:#202124}.title{font-size:21px;font-weight:650}.subtitle{font-size:14px;'
        'fill:#5f6368}.panel{font-size:17px;font-weight:600}.axis{font-size:12px}'
        '.cell{font-size:12px;font-weight:650;fill:#ffffff}.tie{stroke:#202124;'
        'stroke-width:2}.missing{fill:#e8eaed;stroke:#ffffff;stroke-width:1}</style>',
        '<text x="28" y="31" class="title">Pattern-matching algorithm selection</text>',
        '<text x="28" y="53" class="subtitle">Cell color is the lowest median time; dark outline means &lt;3% over runner-up</text>',
    ]

    legend_x = 28
    for algorithm in ALGORITHMS:
        color = COLORS[algorithm]
        svg.append(
            f'<rect x="{legend_x}" y="72" width="20" height="14" fill="{color}"/>'
        )
        svg.append(
            f'<text x="{legend_x+27}" y="84" class="axis">{algorithm}</text>'
        )
        legend_x += 145
    svg.append(
        f'<rect x="{legend_x}" y="72" width="20" height="14" class="missing"/>'
    )
    svg.append(
        f'<text x="{legend_x+27}" y="84" class="axis">unavailable</text>'
    )

    for mode_index, mode in enumerate(modes):
        panel_top = top + mode_index * (panel_height + panel_gap)
        svg.append(
            f'<text x="28" y="{panel_top-16}" class="panel">'
            f'{"First match" if mode == "first" else "All overlapping matches"}</text>'
        )
        for row, k in enumerate(k_values):
            y = panel_top + row * cell_height
            svg.append(
                f'<text x="{left-10}" y="{y+cell_height*0.68:.1f}" '
                f'text-anchor="end" class="axis">{k}</text>'
            )
            for column, patterns in enumerate(pattern_counts):
                x = left + column * cell_width
                winner = winners.get((k, patterns, mode))
                if winner is None:
                    svg.append(
                        f'<rect x="{x}" y="{y}" width="{cell_width}" '
                        f'height="{cell_height}" class="missing"/>'
                    )
                    continue
                algorithm = winner["algorithm"]
                tie_class = (
                    ' class="tie"'
                    if winner["margin"] is not None and winner["margin"] < 3.0
                    else ""
                )
                svg.append(
                    f'<rect x="{x}" y="{y}" width="{cell_width}" '
                    f'height="{cell_height}" fill="{COLORS[algorithm]}" '
                    f'stroke="#ffffff" stroke-width="1"{tie_class}/>'
                )
                svg.append(
                    f'<text x="{x+cell_width/2:.1f}" y="{y+cell_height*0.68:.1f}" '
                    f'text-anchor="middle" class="cell">{LABELS[algorithm]}</text>'
                )

        bottom = panel_top + panel_height
        for column, patterns in enumerate(pattern_counts):
            x = left + (column + 0.5) * cell_width
            svg.append(
                f'<text transform="translate({x:.1f} {bottom+9}) rotate(55)" '
                f'text-anchor="start" class="axis">{patterns:,}</text>'
            )
        svg.append(
            f'<text transform="translate(24 {panel_top+panel_height/2:.1f}) rotate(-90)" '
            'text-anchor="middle" class="axis">pattern length k</text>'
        )
        svg.append(
            f'<text x="{left+len(pattern_counts)*cell_width/2:.1f}" y="{bottom+72}" '
            'text-anchor="middle" class="axis">number of patterns</text>'
        )

    svg.append("</svg>")
    SELECTION_MAP.write_text("\n".join(svg) + "\n")


def main():
    if not RAW.exists() or not STATUS.exists():
        raise SystemExit("Missing algorithm_sweep.csv or cell_status.csv.")
    timings, build_times, validation = load_raw()
    if not timings:
        raise SystemExit("The raw result file contains no successful measurements.")
    statuses = load_statuses()
    validate_algorithms(validation)
    summary = summarize(timings, build_times)
    winners = write_tables(summary)
    write_selection_map(winners, statuses)


if __name__ == "__main__":
    main()
