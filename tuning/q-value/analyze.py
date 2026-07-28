#!/usr/bin/env python3

import csv
import math
import statistics
from collections import defaultdict
from pathlib import Path


ROOT = Path(__file__).resolve().parent
RESULTS = ROOT / "results"
RAW = RESULTS / "q_sweep.csv"
SUMMARY = RESULTS / "q_summary.csv"
BEST = RESULTS / "q_best_by_mode.csv"
PLOT = RESULTS / "q_curves.svg"
OPTIMAL_PLOT = RESULTS / "q_optimal.svg"
MODES = ("first", "all")
COLORS = [
    "#1a73e8",
    "#d93025",
    "#188038",
    "#f29900",
    "#9334e6",
    "#007b83",
    "#e37400",
    "#c5221f",
    "#5f6368",
    "#3f51b5",
    "#00a1c9",
    "#9c6ade",
]


def load_timings():
    timings = defaultdict(dict)
    with RAW.open(newline="") as handle:
        for row in csv.DictReader(handle):
            key = (int(row["k"]), int(row["q"]), row["mode"])
            timings[key][int(row["run"])] = float(row["ns_per_base"])
    return timings


def summarize(timings):
    summary = {}
    for key, values_by_run in timings.items():
        values = list(values_by_run.values())
        summary[key] = {
            "median": statistics.median(values),
            "mean": statistics.fmean(values),
            "stddev": statistics.stdev(values) if len(values) > 1 else 0.0,
            "minimum": min(values),
            "maximum": max(values),
            "runs": len(values),
        }
    return summary


def write_tables(summary):
    best_time = {}
    for (k, _q, mode), values in summary.items():
        best_time[(k, mode)] = min(values["median"], best_time.get((k, mode), float("inf")))

    with SUMMARY.open("w", newline="") as handle:
        writer = csv.writer(handle)
        writer.writerow(
            [
                "k",
                "q",
                "mode",
                "runs",
                "median_ns_per_base",
                "mean_ns_per_base",
                "stddev_ns_per_base",
                "min_ns_per_base",
                "max_ns_per_base",
                "slower_than_best_pct",
            ]
        )
        for (k, q, mode), values in sorted(summary.items()):
            relative = (values["median"] / best_time[(k, mode)] - 1.0) * 100.0
            writer.writerow(
                [
                    k,
                    q,
                    mode,
                    values["runs"],
                    f"{values['median']:.9f}",
                    f"{values['mean']:.9f}",
                    f"{values['stddev']:.9f}",
                    f"{values['minimum']:.9f}",
                    f"{values['maximum']:.9f}",
                    f"{relative:.4f}",
                ]
            )

    rows = []
    for k, mode in sorted(best_time):
        choices = sorted(
            (values["median"], q)
            for (candidate_k, q, candidate_mode), values in summary.items()
            if (candidate_k, candidate_mode) == (k, mode)
        )
        value, q = choices[0]
        second_value, second_q = choices[1] if len(choices) > 1 else choices[0]
        rows.append(
            {
                "k": k,
                "mode": mode,
                "best_q": q,
                "median_ns_per_base": value,
                "second_q": second_q,
                "margin_pct": (second_value / value - 1.0) * 100.0,
            }
        )

    with BEST.open("w", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=[
                "k",
                "mode",
                "best_q",
                "median_ns_per_base",
                "second_q",
                "margin_pct",
            ],
        )
        writer.writeheader()
        for row in rows:
            writer.writerow(
                {
                    **row,
                    "median_ns_per_base": f"{row['median_ns_per_base']:.9f}",
                    "margin_pct": f"{row['margin_pct']:.4f}",
                }
            )

    for legacy_name in (
        "q_best_by_workload.csv",
        "q_recommendation.csv",
        "q_value_tuning.svg",
    ):
        (RESULTS / legacy_name).unlink(missing_ok=True)
    return rows


def polyline(points, color, width=2.0):
    coordinates = " ".join(f"{x:.1f},{y:.1f}" for x, y in points)
    return (
        f'<polyline points="{coordinates}" fill="none" stroke="{color}" '
        f'stroke-width="{width}" stroke-linejoin="round" stroke-linecap="round"/>'
    )


def write_svg(summary):
    width, height = 1280, 950
    left, right = 82, 34
    panel_height = 280
    panel_tops = {"first": 165, "all": 570}
    k_values = sorted({key[0] for key in summary})
    q_values = sorted({key[1] for key in summary})
    k_min, k_max = min(k_values), max(k_values)
    maximum = max(values["median"] for values in summary.values()) * 1.05
    tick_step = maximum / 5.0

    def x(k):
        return left + (k - k_min) / max(1, k_max - k_min) * (
            width - left - right
        )

    def y(value, mode):
        bottom = panel_tops[mode] + panel_height
        return bottom - value / maximum * panel_height

    svg = [
        f'<svg xmlns="http://www.w3.org/2000/svg" width="{width}" height="{height}" '
        f'viewBox="0 0 {width} {height}">',
        '<rect width="100%" height="100%" fill="#ffffff"/>',
        '<style>text{font-family:-apple-system,BlinkMacSystemFont,"Segoe UI",sans-serif;'
        'fill:#202124}.axis{stroke:#5f6368;stroke-width:1}.grid{stroke:#dadce0;'
        'stroke-width:1}.label{font-size:13px}.panel-title{font-size:17px;font-weight:600}'
        '.title{font-size:21px;font-weight:650}.subtitle{font-size:14px;fill:#5f6368}</style>',
        '<text x="82" y="31" class="title">BNDMq q-value tuning</text>',
        '<text x="82" y="53" class="subtitle">Median search time across the mixed corpus and all tested patterns</text>',
    ]

    legend_y = 82
    for index, q in enumerate(q_values):
        column = index % 6
        row = index // 6
        legend_x = left + column * 118
        current_y = legend_y + row * 24
        color = COLORS[q - 1]
        svg.append(
            f'<line x1="{legend_x}" y1="{current_y}" x2="{legend_x+28}" '
            f'y2="{current_y}" stroke="{color}" stroke-width="3"/>'
        )
        svg.append(
            f'<text x="{legend_x+36}" y="{current_y+4}" class="label">q={q}</text>'
        )

    titles = {
        "first": "First-match search",
        "all": "All-match enumeration",
    }
    for mode in MODES:
        top = panel_tops[mode]
        bottom = top + panel_height
        svg.append(f'<text x="{left}" y="{top-17}" class="panel-title">{titles[mode]}</text>')
        for tick in range(6):
            value = tick * tick_step
            current_y = y(value, mode)
            svg.append(
                f'<line x1="{left}" y1="{current_y:.1f}" x2="{width-right}" '
                f'y2="{current_y:.1f}" class="grid"/>'
            )
            svg.append(
                f'<text x="{left-12}" y="{current_y+4:.1f}" text-anchor="end" '
                f'class="label">{value:.2f}</text>'
            )
        svg.append(
            f'<line x1="{left}" y1="{top}" x2="{left}" y2="{bottom}" class="axis"/>'
        )
        svg.append(
            f'<line x1="{left}" y1="{bottom}" x2="{width-right}" y2="{bottom}" class="axis"/>'
        )
        svg.append(
            f'<text transform="translate(25 {(top+bottom)/2}) rotate(-90)" '
            'text-anchor="middle" class="label">median ns/base</text>'
        )

        for q in q_values:
            points = [
                (x(k), y(summary[(k, q, mode)]["median"], mode))
                for k in k_values
                if (k, q, mode) in summary
            ]
            svg.append(polyline(points, COLORS[q - 1]))

        for k in range(((k_min + 9) // 10) * 10, k_max + 1, 10):
            current_x = x(k)
            svg.append(
                f'<line x1="{current_x:.1f}" y1="{bottom}" x2="{current_x:.1f}" '
                f'y2="{bottom+5}" class="axis"/>'
            )
            svg.append(
                f'<text x="{current_x:.1f}" y="{bottom+22}" text-anchor="middle" '
                f'class="label">{k}</text>'
            )
        svg.append(
            f'<text x="{(left+width-right)/2:.1f}" y="{bottom+46}" text-anchor="middle" '
            'class="label">pattern length k</text>'
        )

    svg.append("</svg>")
    PLOT.write_text("\n".join(svg) + "\n")


def write_optimal_svg(rows):
    width, height = 1280, 790
    left, right = 82, 34
    top_1, bottom_1 = 105, 420
    top_2, bottom_2 = 535, 735
    k_values = sorted({row["k"] for row in rows})
    k_min, k_max = min(k_values), max(k_values)
    max_q = 12
    maximum_margin = max(row["margin_pct"] for row in rows)
    margin_axis_min = 0.0
    margin_axis_max = max(5.0, math.ceil(maximum_margin / 5.0) * 5.0)
    mode_colors = {"first": "#1a73e8", "all": "#d93025"}
    mode_labels = {"first": "first match", "all": "all matches"}

    def x(k):
        return left + (k - k_min) / max(1, k_max - k_min) * (
            width - left - right
        )

    def y_q(q):
        return bottom_1 - (q - 1) / (max_q - 1) * (bottom_1 - top_1)

    def y_margin(value):
        return bottom_2 - (value - margin_axis_min) / (
            margin_axis_max - margin_axis_min
        ) * (bottom_2 - top_2)

    svg = [
        f'<svg xmlns="http://www.w3.org/2000/svg" width="{width}" height="{height}" '
        f'viewBox="0 0 {width} {height}">',
        f'<defs><clipPath id="gap-clip"><rect x="{left}" y="{top_2}" '
        f'width="{width-left-right}" height="{bottom_2-top_2}"/></clipPath></defs>',
        '<rect width="100%" height="100%" fill="#ffffff"/>',
        '<style>text{font-family:-apple-system,BlinkMacSystemFont,"Segoe UI",sans-serif;'
        'fill:#202124}.axis{stroke:#5f6368;stroke-width:1}.grid{stroke:#dadce0;'
        'stroke-width:1}.label{font-size:13px}.panel-title{font-size:17px;font-weight:600}'
        '.title{font-size:21px;font-weight:650}.subtitle{font-size:14px;fill:#5f6368}</style>',
        '<text x="82" y="31" class="title">Optimal BNDMq q by pattern length</text>',
        '<text x="82" y="53" class="subtitle">Mixed corpus; the shaded band marks a practical tie with the runner-up</text>',
    ]

    for index, mode in enumerate(MODES):
        legend_x = 790 + index * 170
        color = mode_colors[mode]
        svg.append(
            f'<line x1="{legend_x}" y1="42" x2="{legend_x+30}" y2="42" '
            f'stroke="{color}" stroke-width="3"/>'
        )
        svg.append(
            f'<circle cx="{legend_x+15}" cy="42" r="3.5" fill="{color}"/>'
        )
        svg.append(
            f'<text x="{legend_x+39}" y="46" class="label">{mode_labels[mode]}</text>'
        )

    svg.append(f'<text x="{left}" y="{top_1-17}" class="panel-title">Selected q</text>')
    for q in range(1, max_q + 1):
        current_y = y_q(q)
        svg.append(
            f'<line x1="{left}" y1="{current_y:.1f}" x2="{width-right}" '
            f'y2="{current_y:.1f}" class="grid"/>'
        )
        svg.append(
            f'<text x="{left-12}" y="{current_y+4:.1f}" text-anchor="end" '
            f'class="label">{q}</text>'
        )
    svg.append(
        f'<line x1="{left}" y1="{top_1}" x2="{left}" y2="{bottom_1}" class="axis"/>'
    )
    svg.append(
        f'<line x1="{left}" y1="{bottom_1}" x2="{width-right}" y2="{bottom_1}" class="axis"/>'
    )

    for mode in MODES:
        mode_rows = sorted(
            (row for row in rows if row["mode"] == mode),
            key=lambda row: row["k"],
        )
        points = [(x(row["k"]), y_q(row["best_q"])) for row in mode_rows]
        color = mode_colors[mode]
        svg.append(polyline(points, color, width=2.5))
        for point_x, point_y in points:
            svg.append(
                f'<circle cx="{point_x:.1f}" cy="{point_y:.1f}" r="2.8" '
                f'fill="{color}"/>'
            )

    svg.append(
        f'<text x="{left}" y="{top_2-17}" class="panel-title">'
        "Advantage over second-best q</text>"
    )
    near_tie_top = y_margin(1.0)
    near_tie_bottom = y_margin(0.0)
    svg.append(
        f'<rect x="{left}" y="{near_tie_top:.1f}" width="{width-left-right}" '
        f'height="{near_tie_bottom-near_tie_top:.1f}" fill="#f1f3f4"/>'
    )
    margin_tick = 5.0
    tick = margin_axis_min
    while tick <= margin_axis_max + 1e-9:
        current_y = y_margin(tick)
        svg.append(
            f'<line x1="{left}" y1="{current_y:.1f}" x2="{width-right}" '
            f'y2="{current_y:.1f}" class="grid"/>'
        )
        svg.append(
            f'<text x="{left-12}" y="{current_y+4:.1f}" text-anchor="end" '
            f'class="label">{tick:.0f}%</text>'
        )
        tick += margin_tick
    svg.append(
        f'<line x1="{left}" y1="{top_2}" x2="{left}" y2="{bottom_2}" class="axis"/>'
    )
    svg.append(
        f'<line x1="{left}" y1="{bottom_2}" x2="{width-right}" y2="{bottom_2}" class="axis"/>'
    )

    svg.append('<g clip-path="url(#gap-clip)">')
    for mode in MODES:
        mode_rows = sorted(
            (row for row in rows if row["mode"] == mode),
            key=lambda row: row["k"],
        )
        points = [(x(row["k"]), y_margin(row["margin_pct"])) for row in mode_rows]
        color = mode_colors[mode]
        svg.append(polyline(points, color, width=2.0))
        for point_x, point_y in points:
            svg.append(
                f'<circle cx="{point_x:.1f}" cy="{point_y:.1f}" r="2.5" '
                f'fill="{color}"/>'
            )
    svg.append("</g>")

    for k in range(((k_min + 9) // 10) * 10, k_max + 1, 10):
        current_x = x(k)
        for bottom in (bottom_1, bottom_2):
            svg.append(
                f'<line x1="{current_x:.1f}" y1="{bottom}" x2="{current_x:.1f}" '
                f'y2="{bottom+5}" class="axis"/>'
            )
        svg.append(
            f'<text x="{current_x:.1f}" y="{bottom_2+22}" text-anchor="middle" '
            f'class="label">{k}</text>'
        )
    svg.append(
        f'<text x="{(left+width-right)/2:.1f}" y="{bottom_2+47}" '
        'text-anchor="middle" class="label">pattern length k</text>'
    )
    svg.append("</svg>")
    OPTIMAL_PLOT.write_text("\n".join(svg) + "\n")


def main():
    if not RAW.exists():
        raise SystemExit(f"Missing raw results: {RAW}")
    timings = load_timings()
    if not timings:
        raise SystemExit("Raw result file contains no measurements.")
    summary = summarize(timings)
    best_rows = write_tables(summary)
    write_svg(summary)
    write_optimal_svg(best_rows)


if __name__ == "__main__":
    main()
