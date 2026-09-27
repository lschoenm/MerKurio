#!/usr/bin/env python3

import argparse
import csv
import math
import subprocess
import sys
from collections import defaultdict
from pathlib import Path


ROOT = Path(__file__).resolve().parent
RESULTS = ROOT / "results"
REFINEMENT = RESULTS / "refinement"
PARTS = REFINEMENT / "parts"
BASE_METADATA = RESULTS / "metadata.txt"
BASE_RUN_TIMESTAMP = REFINEMENT / "base_run_timestamp.txt"
WINNERS = RESULTS / "algorithm_winners.csv"
REFINEMENT_RAW = REFINEMENT / "algorithm_sweep.csv"
REFINEMENT_STATUS = REFINEMENT / "cell_status.csv"
MANIFEST = REFINEMENT / "manifest.csv"
ANALYZER = ROOT / "analyze.py"
ALGORITHMS = ("bndmq", "hash", "aho_corasick")
CONFIG_FLAGS = {
    "sequence_length": "--sequence-length",
    "total_bases_per_cell": "--total-bases",
    "target_patterns_per_cell": "--target-patterns-per-cell",
    "max_banks": "--max-banks",
    "runs": "--runs",
    "max_sample_ms": "--max-sample-ms",
    "cell_timeout_seconds": "--cell-timeout-seconds",
    "seed": "--seed",
}


def parse_args():
    parser = argparse.ArgumentParser(
        description=(
            "Refine adjacent algorithm-winner transitions along pattern count "
            "and pattern length with one midpoint measurement wave."
        )
    )
    parser.add_argument("--binary", type=Path, required=True)
    parser.add_argument(
        "--waves",
        type=int,
        default=1,
        help="Adaptive midpoint waves to run [default: 1]",
    )
    parser.add_argument(
        "--runs",
        type=int,
        help="Override the timing rounds recorded in the base metadata",
    )
    parser.add_argument(
        "--dry-run",
        action="store_true",
        help="Print the next refinement cells without running them",
    )
    parser.add_argument(
        "--allow-different-system",
        action="store_true",
        help="Allow refinement on a system different from the base sweep",
    )
    args = parser.parse_args()
    if args.waves < 1:
        parser.error("--waves must be positive")
    if args.runs is not None and args.runs < 1:
        parser.error("--runs must be positive")
    return args


def read_metadata(path=BASE_METADATA):
    metadata = {}
    with path.open() as handle:
        for line in handle:
            key, separator, value = line.rstrip("\n").partition("=")
            if separator:
                metadata[key] = value
    missing = sorted(set(CONFIG_FLAGS) - metadata.keys())
    if missing:
        raise SystemExit(
            "Base metadata is missing refinement settings: " + ", ".join(missing)
        )
    return metadata


def validate_provenance(metadata, allow_different_system):
    current_system = subprocess.run(
        ("uname", "-a"),
        check=True,
        capture_output=True,
        text=True,
    ).stdout.strip()
    if (
        metadata.get("system") != current_system
        and not allow_different_system
    ):
        raise SystemExit(
            "The base sweep was produced on a different system. Run the "
            "refinement there, or use --allow-different-system only if this "
            "difference is intentional."
        )

    timestamp = metadata.get("unix_timestamp")
    if not timestamp:
        raise SystemExit("Base metadata does not contain unix_timestamp.")
    if BASE_RUN_TIMESTAMP.exists():
        recorded = BASE_RUN_TIMESTAMP.read_text().strip()
        if recorded != timestamp:
            raise SystemExit(
                "Existing refinement parts belong to a different base sweep. "
                "Archive or remove results/refinement before starting again."
            )
    else:
        REFINEMENT.mkdir(parents=True, exist_ok=True)
        BASE_RUN_TIMESTAMP.write_text(timestamp + "\n")


def read_winners(path=WINNERS):
    grouped = defaultdict(list)
    with path.open(newline="") as handle:
        for row in csv.DictReader(handle):
            grouped[(int(row["k"]), row["mode"])].append(
                {
                    "patterns": int(row["patterns"]),
                    "algorithm": row["best_algorithm"],
                }
            )
    for rows in grouped.values():
        rows.sort(key=lambda row: row["patterns"])
    return grouped


def geometric_midpoint(lower, upper):
    midpoint = round(math.sqrt(lower * upper))
    return min(upper - 1, max(lower + 1, midpoint))


def transition_targets(grouped_winners):
    targets = []
    for (k, mode), rows in sorted(grouped_winners.items()):
        for left, right in zip(rows, rows[1:]):
            lower = left["patterns"]
            upper = right["patterns"]
            if (
                left["algorithm"] == right["algorithm"]
                or upper - lower <= 1
            ):
                continue
            targets.append(
                {
                    "axis": "patterns",
                    "k": k,
                    "mode": mode,
                    "patterns": geometric_midpoint(lower, upper),
                    "lower_k": k,
                    "upper_k": k,
                    "lower_patterns": lower,
                    "upper_patterns": upper,
                    "lower_algorithm": left["algorithm"],
                    "upper_algorithm": right["algorithm"],
                }
            )

    modes = sorted({mode for _k, mode in grouped_winners})
    for mode in modes:
        k_values = sorted(k for k, candidate_mode in grouped_winners if candidate_mode == mode)
        for lower_k, upper_k in zip(k_values, k_values[1:]):
            if upper_k - lower_k <= 1:
                continue
            midpoint_k = 65 if lower_k <= 64 < upper_k else (lower_k + upper_k) // 2
            lower_by_patterns = {
                row["patterns"]: row["algorithm"]
                for row in grouped_winners[(lower_k, mode)]
            }
            upper_by_patterns = {
                row["patterns"]: row["algorithm"]
                for row in grouped_winners[(upper_k, mode)]
            }
            for patterns in sorted(lower_by_patterns.keys() & upper_by_patterns.keys()):
                lower_algorithm = lower_by_patterns[patterns]
                upper_algorithm = upper_by_patterns[patterns]
                if lower_algorithm == upper_algorithm:
                    continue
                targets.append(
                    {
                        "axis": "k",
                        "k": midpoint_k,
                        "mode": mode,
                        "patterns": patterns,
                        "lower_k": lower_k,
                        "upper_k": upper_k,
                        "lower_patterns": patterns,
                        "upper_patterns": patterns,
                        "lower_algorithm": lower_algorithm,
                        "upper_algorithm": upper_algorithm,
                    }
                )

    unique = {}
    for target in targets:
        key = (target["k"], target["patterns"], target["mode"])
        unique[key] = target
    return [unique[key] for key in sorted(unique)]


def part_directory(target):
    return PARTS / (
        f"{target['mode']}-k{target['k']}-p{target['patterns']}"
    )


def completed_part(target):
    directory = part_directory(target)
    status_path = directory / "cell_status.csv"
    raw_path = directory / "algorithm_sweep.csv"
    if not status_path.exists() or not raw_path.exists():
        return False
    with status_path.open(newline="") as handle:
        algorithms = {row["algorithm"] for row in csv.DictReader(handle)}
    return algorithms == set(ALGORITHMS)


def run_cell(binary, target, metadata, runs):
    directory = part_directory(target)
    directory.mkdir(parents=True, exist_ok=True)
    command = [
        str(binary),
        "--k-values",
        str(target["k"]),
        "--pattern-counts",
        str(target["patterns"]),
        "--algorithms",
        ",".join(ALGORITHMS),
        "--modes",
        target["mode"],
    ]
    for key, flag in CONFIG_FLAGS.items():
        value = str(runs) if key == "runs" else metadata[key]
        command.extend((flag, value))
    command.extend(
        (
            "--output",
            str(directory / "algorithm_sweep.csv"),
            "--status-output",
            str(directory / "cell_status.csv"),
            "--metadata",
            str(directory / "metadata.txt"),
        )
    )
    print(
        f"refining {target['mode']} k={target['k']} "
        f"patterns={target['patterns']} "
        f"({target['axis']} transition: "
        f"k={target['lower_k']}..{target['upper_k']}, "
        f"patterns={target['lower_patterns']}..{target['upper_patterns']})",
        flush=True,
    )
    subprocess.run(command, check=True)


def merge_csv(paths, destination):
    header = None
    rows = []
    for path in sorted(paths):
        with path.open(newline="") as handle:
            reader = csv.reader(handle)
            current_header = next(reader, None)
            if current_header is None:
                continue
            if header is None:
                header = current_header
            elif current_header != header:
                raise SystemExit(f"Incompatible CSV header in {path}")
            rows.extend(reader)
    if header is None:
        return
    destination.parent.mkdir(parents=True, exist_ok=True)
    temporary = destination.with_suffix(destination.suffix + ".tmp")
    with temporary.open("w", newline="") as handle:
        writer = csv.writer(handle, lineterminator="\n")
        writer.writerow(header)
        writer.writerows(rows)
    temporary.replace(destination)


def rebuild_aggregates():
    merge_csv(PARTS.glob("*/algorithm_sweep.csv"), REFINEMENT_RAW)
    merge_csv(PARTS.glob("*/cell_status.csv"), REFINEMENT_STATUS)


def read_manifest():
    if not MANIFEST.exists():
        return {}
    with MANIFEST.open(newline="") as handle:
        return {
            (int(row["k"]), int(row["patterns"]), row["mode"]): row
            for row in csv.DictReader(handle)
        }


def write_manifest(targets, wave):
    entries = read_manifest()
    for target in targets:
        key = (target["k"], target["patterns"], target["mode"])
        entries[key] = {
            **{field: str(value) for field, value in target.items()},
            "wave": str(wave),
            "complete": "yes" if completed_part(target) else "no",
        }
    MANIFEST.parent.mkdir(parents=True, exist_ok=True)
    with MANIFEST.open("w", newline="") as handle:
        fields = [
            "wave",
            "axis",
            "k",
            "patterns",
            "mode",
            "lower_k",
            "upper_k",
            "lower_patterns",
            "upper_patterns",
            "lower_algorithm",
            "upper_algorithm",
            "complete",
        ]
        writer = csv.DictWriter(handle, fieldnames=fields, lineterminator="\n")
        writer.writeheader()
        for key in sorted(entries):
            writer.writerow(entries[key])


def analyze():
    subprocess.run((sys.executable, str(ANALYZER)), check=True)


def main():
    args = parse_args()
    metadata = read_metadata()
    if not args.dry_run:
        if not args.binary.is_file():
            raise SystemExit(f"Benchmark binary does not exist: {args.binary}")
        validate_provenance(metadata, args.allow_different_system)
    runs = args.runs if args.runs is not None else int(metadata["runs"])

    if REFINEMENT_RAW.exists() and REFINEMENT_STATUS.exists():
        analyze()

    for wave in range(1, args.waves + 1):
        targets = transition_targets(read_winners())
        pending = [target for target in targets if not completed_part(target)]
        if not pending:
            if targets and any(completed_part(target) for target in targets):
                write_manifest(targets, wave)
                rebuild_aggregates()
                analyze()
            print("No unmeasured transition midpoints remain.")
            break
        print(
            f"refinement wave {wave}: {len(pending)} cells, "
            f"{len(pending) * len(ALGORITHMS)} algorithm comparisons"
        )
        if args.dry_run:
            for target in pending:
                print(
                    f"{target['mode']},k={target['k']},"
                    f"patterns={target['patterns']},"
                    f"axis={target['axis']},"
                    f"k_bracket={target['lower_k']}..{target['upper_k']},"
                    f"pattern_bracket="
                    f"{target['lower_patterns']}..{target['upper_patterns']}"
                )
            return
        for target in pending:
            run_cell(args.binary, target, metadata, runs)
        write_manifest(targets, wave)
        rebuild_aggregates()
        analyze()

    if not args.dry_run and REFINEMENT_RAW.exists():
        print(f"Refinement results written to {REFINEMENT}")
        print(f"Updated plots written to {RESULTS}")


if __name__ == "__main__":
    main()
