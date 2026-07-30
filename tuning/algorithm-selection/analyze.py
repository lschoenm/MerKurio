#!/usr/bin/env python3

import csv
import math
import random
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
CROSSOVER_PLOT = RESULTS / "crossover_curves.svg"
CORPUS_REGIONS = RESULTS / "corpus_size_regions.csv"
CORPUS_CROSSOVERS = RESULTS / "corpus_size_crossovers.csv"
CORPUS_CROSSOVER_PLOT = RESULTS / "corpus_size_crossovers.svg"
CORPUS_MODELS = RESULTS / "corpus_size_models.csv"
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
    search_timings = defaultdict(list)
    build_times = defaultdict(list)
    corpus_model_samples = defaultdict(list)
    validation = defaultdict(dict)
    reference_bases = set()
    with RAW.open(newline="") as handle:
        for row in csv.DictReader(handle):
            key = (
                int(row["k"]),
                int(row["patterns"]),
                row["mode"],
                row["algorithm"],
            )
            timings[key].append(float(row["reference_ns_per_base"]))
            search_timings[key].append(float(row["search_ns_per_base"]))
            build_times[key].append(float(row["build_ns_per_matcher"]) / 1_000_000.0)
            corpus_model_samples[key].append(
                (
                    float(row["fixed_ns_per_matcher"]),
                    float(row["search_ns_per_base"]),
                )
            )
            reference_bases.add(int(row["target_bases_per_cell"]))
            cell = key[:3]
            observed = (
                int(row["validation_checksum"]),
                int(row["matching_records"]),
                int(row["matches"]),
            )
            previous = validation[cell].setdefault(key[3], observed)
            if previous != observed:
                raise SystemExit(f"Inconsistent validation result within {key}")
    if len(reference_bases) != 1:
        raise SystemExit("Raw results contain inconsistent target corpus sizes.")
    return (
        timings,
        search_timings,
        build_times,
        corpus_model_samples,
        validation,
        reference_bases.pop(),
    )


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


def summarize(timings, search_timings, build_times):
    summary = {}
    for key, values in timings.items():
        summary[key] = {
            "runs": len(values),
            "median": statistics.median(values),
            "mean": statistics.fmean(values),
            "stddev": statistics.stdev(values) if len(values) > 1 else 0.0,
            "minimum": min(values),
            "maximum": max(values),
            "search_median": statistics.median(search_timings[key]),
            "build_ms": statistics.median(build_times[key]),
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
        writer = csv.writer(handle, lineterminator="\n")
        writer.writerow(
            [
                "k",
                "patterns",
                "mode",
                "algorithm",
                "runs",
                "median_build_ms_per_matcher",
                "median_search_ns_per_base",
                "median_reference_ns_per_base",
                "mean_reference_ns_per_base",
                "stddev_reference_ns_per_base",
                "min_reference_ns_per_base",
                "max_reference_ns_per_base",
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
                    f"{values['build_ms']:.6f}",
                    f"{values['search_median']:.9f}",
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
        writer = csv.writer(handle, lineterminator="\n")
        writer.writerow(
            [
                "k",
                "patterns",
                "mode",
                "best_algorithm",
                "median_reference_ns_per_base",
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


def write_selection_map(winners, statuses, reference_bases):
    k_values = sorted({key[0] for key in statuses})
    pattern_counts = sorted({key[1] for key in statuses})
    modes = [mode for mode in MODES if any(key[2] == mode for key in statuses)]
    cell_width = max(58, min(92, 1120 // max(1, len(pattern_counts))))
    cell_height = max(25, min(38, 500 // max(1, len(k_values))))
    left = 96
    right = 28
    top = 130
    panel_gap = 130
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
        '<text x="28" y="53" class="subtitle">Cell color is the lowest median '
        f'{compact_bases(reference_bases)}-base reference time; dark outline '
        'means &lt;3% over runner-up</text>',
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


def polyline(points, color, width=2.0):
    coordinates = " ".join(f"{x:.1f},{y:.1f}" for x, y in points)
    return (
        f'<polyline points="{coordinates}" fill="none" stroke="{color}" '
        f'stroke-width="{width}" stroke-linejoin="round" stroke-linecap="round"/>'
    )


def write_crossover_plot(summary, statuses, reference_bases):
    k_values = sorted({key[0] for key in statuses})
    pattern_counts = sorted({key[1] for key in statuses})
    modes = [mode for mode in MODES if any(key[2] == mode for key in statuses)]
    columns = min(3, max(1, len(k_values)))
    rows_per_mode = math.ceil(len(k_values) / columns)
    panel_width, panel_height = 340, 205
    left, right, top = 66, 28, 125
    mode_heading = 38
    mode_gap = 45
    width = max(760, left + columns * panel_width + right)
    height = (
        top
        + len(modes) * (mode_heading + rows_per_mode * panel_height)
        + max(0, len(modes) - 1) * mode_gap
        + 35
    )

    relative = {}
    for k in k_values:
        for patterns in pattern_counts:
            for mode in modes:
                available = {
                    algorithm: summary[(k, patterns, mode, algorithm)]["median"]
                    for algorithm in ALGORITHMS
                    if (k, patterns, mode, algorithm) in summary
                }
                if available:
                    best = min(available.values())
                    for algorithm, value in available.items():
                        relative[(k, patterns, mode, algorithm)] = value / best

    log_min = math.log10(min(pattern_counts))
    log_max = math.log10(max(pattern_counts))
    global_max = max(relative.values(), default=1.0)
    y_max = max(1.25, min(5.0, math.ceil(global_max * 2.0) / 2.0))
    clipping_note = (
        f"; values above {y_max:.1f}× are clipped"
        if global_max > y_max
        else ""
    )

    svg = [
        f'<svg xmlns="http://www.w3.org/2000/svg" width="{width}" height="{height}" '
        f'viewBox="0 0 {width} {height}">',
        '<rect width="100%" height="100%" fill="#ffffff"/>',
        '<style>text{font-family:-apple-system,BlinkMacSystemFont,"Segoe UI",sans-serif;'
        'fill:#202124}.title{font-size:21px;font-weight:650}.subtitle{font-size:14px;'
        'fill:#5f6368}.mode{font-size:18px;font-weight:600}.panel{font-size:14px;'
        'font-weight:600}.axis{font-size:11px}.grid{stroke:#dadce0;stroke-width:1}'
        '.frame{stroke:#5f6368;stroke-width:1;fill:none}</style>',
        '<text x="28" y="31" class="title">Algorithm crossover curves</text>',
        '<text x="28" y="53" class="subtitle">Reference runtime at '
        f'{compact_bases(reference_bases)} bases relative to the fastest available algorithm'
        f"{clipping_note}</text>",
    ]
    legend_x = 28
    for algorithm in ALGORITHMS:
        svg.append(
            f'<line x1="{legend_x}" y1="79" x2="{legend_x+27}" y2="79" '
            f'stroke="{COLORS[algorithm]}" stroke-width="3"/>'
        )
        svg.append(
            f'<text x="{legend_x+34}" y="83" class="axis">{algorithm}</text>'
        )
        legend_x += 155

    for mode_index, mode in enumerate(modes):
        mode_top = top + mode_index * (
            mode_heading + rows_per_mode * panel_height + mode_gap
        )
        svg.append(
            f'<text x="28" y="{mode_top}" class="mode">'
            f'{"First match" if mode == "first" else "All overlapping matches"}</text>'
        )
        for index, k in enumerate(k_values):
            row, column = divmod(index, columns)
            panel_x = left + column * panel_width
            panel_y = mode_top + mode_heading + row * panel_height
            plot_left = panel_x + 42
            plot_right = panel_x + panel_width - 18
            plot_top = panel_y + 24
            plot_bottom = panel_y + panel_height - 42

            def x(patterns):
                if log_max == log_min:
                    return (plot_left + plot_right) / 2
                return plot_left + (
                    (math.log10(patterns) - log_min) / (log_max - log_min)
                ) * (plot_right - plot_left)

            def y(value):
                clipped = min(value, y_max)
                return plot_bottom - (clipped - 1.0) / (y_max - 1.0) * (
                    plot_bottom - plot_top
                )

            svg.append(
                f'<text x="{plot_left}" y="{panel_y+15}" class="panel">k={k}</text>'
            )
            for tick in (1.0, (1.0 + y_max) / 2.0, y_max):
                current_y = y(tick)
                svg.append(
                    f'<line x1="{plot_left}" y1="{current_y:.1f}" '
                    f'x2="{plot_right}" y2="{current_y:.1f}" class="grid"/>'
                )
                svg.append(
                    f'<text x="{plot_left-7}" y="{current_y+4:.1f}" '
                    f'text-anchor="end" class="axis">{tick:.1f}×</text>'
                )
            svg.append(
                f'<rect x="{plot_left}" y="{plot_top}" '
                f'width="{plot_right-plot_left}" height="{plot_bottom-plot_top}" '
                'class="frame"/>'
            )

            for algorithm in ALGORITHMS:
                segments = []
                current = []
                for patterns in pattern_counts:
                    key = (k, patterns, mode, algorithm)
                    if key in relative:
                        current.append((x(patterns), y(relative[key])))
                    elif current:
                        segments.append(current)
                        current = []
                if current:
                    segments.append(current)
                for points in segments:
                    if len(points) > 1:
                        svg.append(polyline(points, COLORS[algorithm], 2.0))
                    for point_x, point_y in points:
                        svg.append(
                            f'<circle cx="{point_x:.1f}" cy="{point_y:.1f}" '
                            f'r="2.3" fill="{COLORS[algorithm]}"/>'
                        )

            tick_counts = sorted(
                {
                    pattern_counts[0],
                    pattern_counts[len(pattern_counts) // 2],
                    pattern_counts[-1],
                }
            )
            for patterns in tick_counts:
                current_x = x(patterns)
                svg.append(
                    f'<text transform="translate({current_x:.1f} {plot_bottom+10}) '
                    f'rotate(45)" text-anchor="start" class="axis">{patterns:,}</text>'
                )
            if row == rows_per_mode - 1:
                svg.append(
                    f'<text x="{(plot_left+plot_right)/2:.1f}" '
                    f'y="{panel_y+panel_height-5}" text-anchor="middle" '
                    'class="axis">patterns</text>'
                )

    svg.append("</svg>")
    CROSSOVER_PLOT.write_text("\n".join(svg) + "\n")


def percentile(values, percentile_value):
    if not values:
        return None
    ordered = sorted(values)
    position = (len(ordered) - 1) * percentile_value
    lower = math.floor(position)
    upper = math.ceil(position)
    if lower == upper:
        return ordered[lower]
    fraction = position - lower
    return ordered[lower] * (1.0 - fraction) + ordered[upper] * fraction


def lower_envelope(lines):
    if not lines:
        return []
    boundaries = {0.0}
    choices = list(lines.items())
    for index, (_left_algorithm, (left_fixed, left_rate)) in enumerate(choices):
        for _right_algorithm, (right_fixed, right_rate) in choices[index + 1 :]:
            rate_difference = left_rate - right_rate
            if abs(rate_difference) < 1e-15:
                continue
            crossing = (right_fixed - left_fixed) / rate_difference
            if crossing > 0 and math.isfinite(crossing):
                boundaries.add(crossing)

    ordered = sorted(boundaries)
    interval_winners = []
    for index, lower in enumerate(ordered):
        if index + 1 < len(ordered):
            probe = (lower + ordered[index + 1]) / 2.0
        else:
            probe = max(1.0, lower * 2.0 + 1.0)
        winner = min(
            lines,
            key=lambda algorithm: lines[algorithm][0]
            + lines[algorithm][1] * probe,
        )
        if not interval_winners or interval_winners[-1] != winner:
            interval_winners.append(winner)

    regions = []
    start = 0.0
    for index, algorithm in enumerate(interval_winners):
        if index + 1 < len(interval_winners):
            next_algorithm = interval_winners[index + 1]
            fixed, rate = lines[algorithm]
            next_fixed, next_rate = lines[next_algorithm]
            end = (next_fixed - fixed) / (rate - next_rate)
        else:
            end = None
        regions.append(
            {
                "algorithm": algorithm,
                "start": start,
                "end": end,
            }
        )
        if end is not None:
            start = end
    return regions


def bootstrap_transition_intervals(samples, point_regions, seed, rounds=1000):
    transitions = [
        (region["algorithm"], point_regions[index + 1]["algorithm"])
        for index, region in enumerate(point_regions[:-1])
    ]
    observed = {transition: [] for transition in transitions}
    rng = random.Random(seed)

    for _ in range(rounds):
        lines = {}
        for algorithm, values in samples.items():
            resampled = [values[rng.randrange(len(values))] for _ in values]
            lines[algorithm] = (
                statistics.median(value[0] for value in resampled),
                statistics.median(value[1] for value in resampled),
            )
        regions = lower_envelope(lines)
        for index, region in enumerate(regions[:-1]):
            transition = (region["algorithm"], regions[index + 1]["algorithm"])
            if transition in observed:
                observed[transition].append(region["end"])

    intervals = {}
    for transition, values in observed.items():
        intervals[transition] = {
            "low": percentile(values, 0.025),
            "high": percentile(values, 0.975),
            "support": len(values) / rounds * 100.0,
        }
    return intervals


def compact_bases(value):
    suffixes = ((1e12, "T"), (1e9, "G"), (1e6, "M"), (1e3, "k"))
    for scale, suffix in suffixes:
        if value >= scale:
            scaled = value / scale
            precision = 0 if scaled >= 10 else 1
            return f"{scaled:.{precision}f}{suffix}"
    return f"{value:.0f}"


def corpus_crossover_color(value):
    logarithm = max(2.0, min(12.0, math.log10(max(1.0, value))))
    fraction = (logarithm - 2.0) / 10.0
    start = (254, 229, 153)
    end = (84, 39, 143)
    color = tuple(
        round(left + (right - left) * fraction)
        for left, right in zip(start, end)
    )
    return f"#{color[0]:02x}{color[1]:02x}{color[2]:02x}", fraction


def write_corpus_crossover_plot(regions_by_cell, availability, statuses):
    k_values = sorted({key[0] for key in statuses})
    pattern_counts = sorted({key[1] for key in statuses})
    modes = [mode for mode in MODES if any(key[2] == mode for key in statuses)]
    cell_width = max(82, min(106, 1400 // max(1, len(pattern_counts))))
    cell_height = max(28, min(38, 600 // max(1, len(k_values))))
    left, right, top = 96, 28, 145
    panel_gap = 125
    panel_height = len(k_values) * cell_height
    width = max(880, left + len(pattern_counts) * cell_width + right)
    height = (
        top
        + len(modes) * panel_height
        + max(0, len(modes) - 1) * panel_gap
        + 105
    )

    svg = [
        f'<svg xmlns="http://www.w3.org/2000/svg" width="{width}" height="{height}" '
        f'viewBox="0 0 {width} {height}">',
        '<rect width="100%" height="100%" fill="#ffffff"/>',
        '<style>text{font-family:-apple-system,BlinkMacSystemFont,"Segoe UI",sans-serif;'
        'fill:#202124}.title{font-size:21px;font-weight:650}.subtitle{font-size:14px;'
        'fill:#5f6368}.panel{font-size:17px;font-weight:600}.axis{font-size:11px}'
        '.cell{font-size:10px;font-weight:650}.missing{fill:#e8eaed;stroke:#ffffff;'
        'stroke-width:1}</style>',
        '<text x="28" y="31" class="title">Corpus-size algorithm crossovers (in searched bases)</text>',
        '<text x="28" y="53" class="subtitle">First lower-envelope transition among successful algorithms; +N marks additional transitions</text>',
        '<text x="28" y="82" class="axis">earlier crossover</text>',
    ]
    legend_x = 125
    for index, (value, label) in enumerate(
        ((1e2, "100"), (1e4, "10k"), (1e6, "1M"), (1e8, "100M"), (1e10, "10G"))
    ):
        color, _fraction = corpus_crossover_color(value)
        x = legend_x + index * 75
        svg.append(f'<rect x="{x}" y="69" width="42" height="14" fill="{color}"/>')
        svg.append(f'<text x="{x+21}" y="96" text-anchor="middle" class="axis">{label}</text>')
    svg.append(
        f'<text x="{legend_x+5*75+4}" y="82" class="axis">later crossover · gray = no observed transition</text>'
    )

    for mode_index, mode in enumerate(modes):
        panel_top = top + mode_index * (panel_height + panel_gap)
        svg.append(
            f'<text x="28" y="{panel_top-20}" class="panel">'
            f'{"First match" if mode == "first" else "All overlapping matches"}</text>'
        )
        for row, k in enumerate(k_values):
            y = panel_top + row * cell_height
            svg.append(
                f'<text x="{left-10}" y="{y+cell_height*0.66:.1f}" '
                f'text-anchor="end" class="axis">{k}</text>'
            )
            for column, patterns in enumerate(pattern_counts):
                x = left + column * cell_width
                regions = regions_by_cell.get((k, patterns, mode))
                if not regions:
                    svg.append(
                        f'<rect x="{x}" y="{y}" width="{cell_width}" '
                        f'height="{cell_height}" class="missing"/>'
                    )
                    continue
                if len(regions) == 1:
                    color = "#f1f3f4"
                    successful, eligible = availability[(k, patterns, mode)]
                    if successful == eligible:
                        label = f"{LABELS[regions[0]['algorithm']]} throughout"
                    else:
                        label = (
                            f"{LABELS[regions[0]['algorithm']]} throughout "
                            f"{successful}/{eligible}"
                        )
                    text_color = "#202124"
                else:
                    first = regions[0]
                    second = regions[1]
                    color, fraction = corpus_crossover_color(first["end"])
                    extra = f" +{len(regions)-2}" if len(regions) > 2 else ""
                    label = (
                        f"{LABELS[first['algorithm']]}→"
                        f"{LABELS[second['algorithm']]} {compact_bases(first['end'])}{extra}"
                    )
                    text_color = "#ffffff" if fraction > 0.48 else "#202124"
                svg.append(
                    f'<rect x="{x}" y="{y}" width="{cell_width}" '
                    f'height="{cell_height}" fill="{color}" '
                    'stroke="#ffffff" stroke-width="1"/>'
                )
                svg.append(
                    f'<text x="{x+cell_width/2:.1f}" y="{y+cell_height*0.66:.1f}" '
                    f'text-anchor="middle" class="cell" style="fill:{text_color}">{label}</text>'
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
            f'<text x="{left+len(pattern_counts)*cell_width/2:.1f}" y="{bottom+76}" '
            'text-anchor="middle" class="axis">number of patterns</text>'
        )

    svg.append("</svg>")
    CORPUS_CROSSOVER_PLOT.write_text("\n".join(svg) + "\n")


def write_corpus_crossover_outputs(corpus_model_samples, statuses):
    cells = sorted({key[:3] for key in corpus_model_samples})
    regions_by_cell = {}
    crossover_rows = []
    availability = {}

    with CORPUS_MODELS.open("w", newline="") as handle:
        writer = csv.writer(handle, lineterminator="\n")
        writer.writerow(
            [
                "k",
                "patterns",
                "mode",
                "algorithm",
                "runs",
                "median_fixed_ns_per_matcher",
                "fixed_p2_5_ns_per_matcher",
                "fixed_p97_5_ns_per_matcher",
                "median_search_ns_per_base",
                "search_p2_5_ns_per_base",
                "search_p97_5_ns_per_base",
            ]
        )
        for key, values in sorted(corpus_model_samples.items()):
            fixed = [value[0] for value in values]
            rates = [value[1] for value in values]
            writer.writerow(
                [
                    *key,
                    len(values),
                    f"{statistics.median(fixed):.9f}",
                    f"{percentile(fixed, 0.025):.9f}",
                    f"{percentile(fixed, 0.975):.9f}",
                    f"{statistics.median(rates):.9f}",
                    f"{percentile(rates, 0.025):.9f}",
                    f"{percentile(rates, 0.975):.9f}",
                ]
            )

    with CORPUS_REGIONS.open("w", newline="") as handle:
        writer = csv.writer(handle, lineterminator="\n")
        writer.writerow(
            [
                "k",
                "patterns",
                "mode",
                "region",
                "algorithm",
                "min_bases_inclusive",
                "max_bases_exclusive",
            ]
        )
        for k, patterns, mode in cells:
            samples = {
                algorithm: corpus_model_samples[(k, patterns, mode, algorithm)]
                for algorithm in ALGORITHMS
                if (k, patterns, mode, algorithm) in corpus_model_samples
            }
            lines = {
                algorithm: (
                    statistics.median(value[0] for value in values),
                    statistics.median(value[1] for value in values),
                )
                for algorithm, values in samples.items()
            }
            regions = lower_envelope(lines)
            regions_by_cell[(k, patterns, mode)] = regions
            eligible = sum(
                statuses.get((k, patterns, mode, algorithm), {}).get("status")
                != "invalid"
                for algorithm in ALGORITHMS
            )
            availability[(k, patterns, mode)] = (len(samples), eligible)
            seed = (
                0x6A09_E667
                ^ k * 0x9E37
                ^ patterns * 0x85EB
                ^ (0 if mode == "first" else 0xC2B2_AE35)
            )
            intervals = bootstrap_transition_intervals(samples, regions, seed)
            for region_index, region in enumerate(regions):
                writer.writerow(
                    [
                        k,
                        patterns,
                        mode,
                        region_index,
                        region["algorithm"],
                        f"{region['start']:.6f}",
                        "" if region["end"] is None else f"{region['end']:.6f}",
                    ]
                )
                if region["end"] is not None:
                    next_algorithm = regions[region_index + 1]["algorithm"]
                    interval = intervals[(region["algorithm"], next_algorithm)]
                    crossover_rows.append(
                        [
                            k,
                            patterns,
                            mode,
                            region_index,
                            region["algorithm"],
                            next_algorithm,
                            region["end"],
                            interval["low"],
                            interval["high"],
                            interval["support"],
                        ]
                    )

    with CORPUS_CROSSOVERS.open("w", newline="") as handle:
        writer = csv.writer(handle, lineterminator="\n")
        writer.writerow(
            [
                "k",
                "patterns",
                "mode",
                "transition",
                "from_algorithm",
                "to_algorithm",
                "crossover_bases",
                "bootstrap_ci_low_bases",
                "bootstrap_ci_high_bases",
                "bootstrap_transition_support_pct",
            ]
        )
        for row in crossover_rows:
            writer.writerow(
                [
                    *row[:6],
                    f"{row[6]:.6f}",
                    "" if row[7] is None else f"{row[7]:.6f}",
                    "" if row[8] is None else f"{row[8]:.6f}",
                    f"{row[9]:.2f}",
                ]
            )

    write_corpus_crossover_plot(regions_by_cell, availability, statuses)
    return regions_by_cell, crossover_rows


def main():
    if not RAW.exists() or not STATUS.exists():
        raise SystemExit("Missing algorithm_sweep.csv or cell_status.csv.")
    (
        timings,
        search_timings,
        build_times,
        corpus_model_samples,
        validation,
        reference_bases,
    ) = load_raw()
    if not timings:
        raise SystemExit("The raw result file contains no successful measurements.")
    statuses = load_statuses()
    validate_algorithms(validation)
    summary = summarize(timings, search_timings, build_times)
    winners = write_tables(summary)
    write_selection_map(winners, statuses, reference_bases)
    write_crossover_plot(summary, statuses, reference_bases)
    write_corpus_crossover_outputs(corpus_model_samples, statuses)


if __name__ == "__main__":
    main()
