#!/usr/bin/env python3

import csv
import math
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
REGRET_TABLE = RESULTS / "selector_regret.csv"
REGRET_PLOT = RESULTS / "selector_regret.svg"
SELECTOR_RULES = RESULTS / "selector_rules.txt"
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


def polyline(points, color, width=2.0):
    coordinates = " ".join(f"{x:.1f},{y:.1f}" for x, y in points)
    return (
        f'<polyline points="{coordinates}" fill="none" stroke="{color}" '
        f'stroke-width="{width}" stroke-linejoin="round" stroke-linecap="round"/>'
    )


def write_crossover_plot(summary, statuses):
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
        '<text x="28" y="53" class="subtitle">Runtime relative to the fastest available algorithm in each cell'
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


def split_timing_rounds(timings):
    training = {}
    evaluation = {}
    for key, values in timings.items():
        split = max(1, len(values) // 2)
        training[key] = statistics.median(values[:split])
        held_out = values[split:]
        evaluation[key] = statistics.median(held_out if held_out else values[-1:])
    return training, evaluation


def cell_times_by_mode(timings, mode):
    cells = defaultdict(dict)
    for (k, patterns, candidate_mode, algorithm), value in timings.items():
        if candidate_mode == mode:
            cells[(k, patterns)][algorithm] = value
    return cells


def leaf_choice(samples):
    best_algorithm = None
    best_loss = float("inf")
    for algorithm in ALGORITHMS:
        loss = 0.0
        for _k, _patterns, times in samples:
            oracle = min(times.values())
            selected = times.get(algorithm, oracle * 10.0)
            loss += math.log(max(1.0, selected / oracle))
        if loss < best_loss:
            best_algorithm = algorithm
            best_loss = loss
    return best_algorithm, best_loss


def fit_selector(samples, depth=0, max_depth=4):
    algorithm, base_loss = leaf_choice(samples)
    if depth >= max_depth or len(samples) < 4:
        return {"algorithm": algorithm}

    best_split = None
    for feature_index, feature in ((0, "k"), (1, "patterns")):
        values = sorted({sample[feature_index] for sample in samples})
        for lower, upper in zip(values, values[1:]):
            threshold = (lower + upper) / 2.0
            left = [sample for sample in samples if sample[feature_index] <= threshold]
            right = [sample for sample in samples if sample[feature_index] > threshold]
            if not left or not right:
                continue
            _left_algorithm, left_loss = leaf_choice(left)
            _right_algorithm, right_loss = leaf_choice(right)
            split_loss = left_loss + right_loss
            if best_split is None or split_loss < best_split["loss"]:
                best_split = {
                    "feature": feature,
                    "threshold": threshold,
                    "left": left,
                    "right": right,
                    "loss": split_loss,
                }

    if best_split is None or best_split["loss"] >= base_loss - 1e-9:
        return {"algorithm": algorithm}
    return {
        "feature": best_split["feature"],
        "threshold": best_split["threshold"],
        "left": fit_selector(best_split["left"], depth + 1, max_depth),
        "right": fit_selector(best_split["right"], depth + 1, max_depth),
    }


def predict_selector(tree, k, patterns):
    while "algorithm" not in tree:
        value = k if tree["feature"] == "k" else patterns
        tree = tree["left"] if value <= tree["threshold"] else tree["right"]
    return tree["algorithm"]


def format_selector(tree, indent=""):
    if "algorithm" in tree:
        return [f"{indent}use {tree['algorithm']}"]
    threshold = tree["threshold"]
    threshold_text = (
        str(int(threshold)) if threshold.is_integer() else f"{threshold:.1f}"
    )
    lines = [f"{indent}if {tree['feature']} <= {threshold_text}:"]
    lines.extend(format_selector(tree["left"], indent + "  "))
    lines.append(f"{indent}else:")
    lines.extend(format_selector(tree["right"], indent + "  "))
    return lines


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


def evaluate_selectors(timings, statuses):
    training, evaluation = split_timing_rounds(timings)
    modes = sorted(
        {key[2] for key in statuses},
        key=lambda mode: MODES.index(mode) if mode in MODES else len(MODES),
    )
    trees = {}
    regrets = {}
    rule_lines = [
        "Shallow selector fitted on the first half of timing rounds.",
        "Regret is evaluated on the held-out second half.",
        "Missing/timed-out algorithms receive a 10x training penalty.",
        "",
    ]

    for mode in modes:
        training_cells = cell_times_by_mode(training, mode)
        samples = [
            (k, patterns, times)
            for (k, patterns), times in sorted(training_cells.items())
        ]
        if not samples:
            continue
        tree = fit_selector(samples)
        trees[mode] = tree
        rule_lines.append(f"[{mode}]")
        rule_lines.extend(format_selector(tree))
        rule_lines.append("")

        evaluation_cells = cell_times_by_mode(evaluation, mode)
        for (k, patterns), times in evaluation_cells.items():
            selected = predict_selector(tree, k, patterns)
            oracle_algorithm, oracle_time = min(
                times.items(), key=lambda item: item[1]
            )
            selected_time = times.get(selected)
            regrets[(k, patterns, mode)] = {
                "selected": selected,
                "oracle": oracle_algorithm,
                "regret": (
                    None
                    if selected_time is None
                    else (selected_time / oracle_time - 1.0) * 100.0
                ),
                "selected_time": selected_time,
                "oracle_time": oracle_time,
            }
    SELECTOR_RULES.write_text("\n".join(rule_lines))
    return regrets, trees


def regret_color(value):
    if value is None:
        return "#e8eaed"
    stops = [
        (0.0, (24, 128, 56)),
        (3.0, (251, 188, 4)),
        (10.0, (242, 153, 0)),
        (25.0, (217, 48, 37)),
    ]
    clipped = max(0.0, min(25.0, value))
    for (left_value, left_color), (right_value, right_color) in zip(
        stops, stops[1:]
    ):
        if clipped <= right_value:
            fraction = (clipped - left_value) / (right_value - left_value)
            color = tuple(
                round(left + (right - left) * fraction)
                for left, right in zip(left_color, right_color)
            )
            return f"#{color[0]:02x}{color[1]:02x}{color[2]:02x}"
    return "#d93025"


def write_regret_outputs(regrets, statuses):
    with REGRET_TABLE.open("w", newline="") as handle:
        writer = csv.writer(handle)
        writer.writerow(
            [
                "k",
                "patterns",
                "mode",
                "selected_algorithm",
                "oracle_algorithm",
                "held_out_regret_pct",
                "selected_ns_per_base",
                "oracle_ns_per_base",
            ]
        )
        for (k, patterns, mode), values in sorted(regrets.items()):
            writer.writerow(
                [
                    k,
                    patterns,
                    mode,
                    values["selected"],
                    values["oracle"],
                    (
                        ""
                        if values["regret"] is None
                        else f"{values['regret']:.4f}"
                    ),
                    (
                        ""
                        if values["selected_time"] is None
                        else f"{values['selected_time']:.9f}"
                    ),
                    f"{values['oracle_time']:.9f}",
                ]
            )

    k_values = sorted({key[0] for key in statuses})
    pattern_counts = sorted({key[1] for key in statuses})
    modes = [mode for mode in MODES if any(key[2] == mode for key in statuses)]
    cell_width = max(64, min(94, 1160 // max(1, len(pattern_counts))))
    cell_height = max(27, min(42, 540 // max(1, len(k_values))))
    left, right, top = 96, 28, 145
    panel_gap = 130
    panel_height = len(k_values) * cell_height
    width = max(760, left + len(pattern_counts) * cell_width + right)
    height = (
        top
        + len(modes) * panel_height
        + max(0, len(modes) - 1) * panel_gap
        + 105
    )

    summaries = {}
    for mode in modes:
        values = [
            result["regret"]
            for key, result in regrets.items()
            if key[2] == mode and result["regret"] is not None
        ]
        summaries[mode] = {
            "median": statistics.median(values) if values else None,
            "p95": percentile(values, 0.95),
            "maximum": max(values) if values else None,
        }

    svg = [
        f'<svg xmlns="http://www.w3.org/2000/svg" width="{width}" height="{height}" '
        f'viewBox="0 0 {width} {height}">',
        '<rect width="100%" height="100%" fill="#ffffff"/>',
        '<style>text{font-family:-apple-system,BlinkMacSystemFont,"Segoe UI",sans-serif;'
        'fill:#202124}.title{font-size:21px;font-weight:650}.subtitle{font-size:14px;'
        'fill:#5f6368}.panel{font-size:17px;font-weight:600}.axis{font-size:12px}'
        '.cell{font-size:11px;font-weight:650}.missing{fill:#e8eaed;stroke:#ffffff;'
        'stroke-width:1}</style>',
        '<text x="28" y="31" class="title">Held-out selector regret</text>',
        '<text x="28" y="53" class="subtitle">Shallow decision tree trained on the first half of rounds; slowdown versus held-out cell oracle</text>',
    ]
    legend = [(0, "#188038"), (3, "#fbbc04"), (10, "#f29900"), (25, "#d93025")]
    legend_x = 28
    for value, color in legend:
        svg.append(
            f'<rect x="{legend_x}" y="72" width="28" height="14" fill="{color}"/>'
        )
        svg.append(
            f'<text x="{legend_x+35}" y="84" class="axis">{value}%</text>'
        )
        legend_x += 90
    svg.append(
        f'<rect x="{legend_x}" y="72" width="28" height="14" class="missing"/>'
    )
    svg.append(
        f'<text x="{legend_x+35}" y="84" class="axis">selected algorithm unavailable</text>'
    )

    for mode_index, mode in enumerate(modes):
        panel_top = top + mode_index * (panel_height + panel_gap)
        summary = summaries[mode]
        summary_text = (
            "no evaluable cells"
            if summary["median"] is None
            else (
                f"median {summary['median']:.2f}% · "
                f"p95 {summary['p95']:.2f}% · max {summary['maximum']:.2f}%"
            )
        )
        svg.append(
            f'<text x="28" y="{panel_top-27}" class="panel">'
            f'{"First match" if mode == "first" else "All overlapping matches"}</text>'
        )
        svg.append(
            f'<text x="242" y="{panel_top-27}" class="subtitle">{summary_text}</text>'
        )
        for row, k in enumerate(k_values):
            y = panel_top + row * cell_height
            svg.append(
                f'<text x="{left-10}" y="{y+cell_height*0.66:.1f}" '
                f'text-anchor="end" class="axis">{k}</text>'
            )
            for column, patterns in enumerate(pattern_counts):
                x = left + column * cell_width
                result = regrets.get((k, patterns, mode))
                if result is None:
                    svg.append(
                        f'<rect x="{x}" y="{y}" width="{cell_width}" '
                        f'height="{cell_height}" class="missing"/>'
                    )
                    continue
                value = result["regret"]
                color = regret_color(value)
                svg.append(
                    f'<rect x="{x}" y="{y}" width="{cell_width}" '
                    f'height="{cell_height}" fill="{color}" '
                    'stroke="#ffffff" stroke-width="1"/>'
                )
                label = (
                    f"{LABELS[result['selected']]} ×"
                    if value is None
                    else f"{LABELS[result['selected']]} {value:.1f}"
                )
                text_color = "#ffffff" if value is not None and value >= 10 else "#202124"
                svg.append(
                    f'<text x="{x+cell_width/2:.1f}" y="{y+cell_height*0.66:.1f}" '
                    f'text-anchor="middle" class="cell" style="fill:{text_color}">{label}</text>'
                )

        bottom = panel_top + panel_height
        for column, patterns in enumerate(pattern_counts):
            x = left + (column + 0.5) * cell_width
            svg.append(
                f'<text transform="translate({x:.1f} {bottom+10}) rotate(55)" '
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
    REGRET_PLOT.write_text("\n".join(svg) + "\n")


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
    write_crossover_plot(summary, statuses)
    regrets, _trees = evaluate_selectors(timings, statuses)
    write_regret_outputs(regrets, statuses)


if __name__ == "__main__":
    main()
