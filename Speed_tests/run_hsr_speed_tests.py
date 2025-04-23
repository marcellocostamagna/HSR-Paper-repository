# Script to launche the hsr_speed_test.py script multiple times
# and aggregate the results.

from __future__ import annotations
import argparse, csv, io, os, re, statistics as stats, subprocess, sys
from pathlib import Path
from typing import Dict, List

# ── robust float pattern (allows optional +/‑ sign) ────────────────── )
FLOAT = r"[-+]?\d*\.\d+(?:[eE][-+]?\d+)?|[-+]?\d+(?:\.\d*)?"

PATTERNS: Dict[str, re.Pattern[str]] = {
    "molecules"   : re.compile(r"Molecules processed\s*:\s*([\d,]+)"),
    "avg_fp"      : re.compile(rf"Average fingerprint time\s*:\s*({FLOAT})"),
    "max_fp"      : re.compile(rf"Maximum fingerprint time\s*:\s*({FLOAT})"),
    "min_fp"      : re.compile(rf"Minimum fingerprint time\s*:\s*({FLOAT})"),
    "seq_fp_time" : re.compile(rf"Sequential fp time\s*:\s*({FLOAT})"),
    "sim_count"   : re.compile(r"([0-9,]+)\s+similarity comparisons", re.I),
    "sim_time"    : re.compile(rf"similarity comparisons.*in\s*({FLOAT})\s*s", re.I),
    "wall_fp"     : re.compile(rf"Wall‑clock fp phase\s*:\s*({FLOAT})"),
    "wall_dl"     : re.compile(rf"Wall‑clock download phase\s*:\s*({FLOAT})"),
    "runtime"     : re.compile(rf"End‑to‑end runtime\s*:\s*({FLOAT})"),
}

# ╭──────────────────────────────────────────────────────────╮
# │ Extract metrics from a single run’s stdout               │
# ╰──────────────────────────────────────────────────────────╯
def parse_stdout(text: str) -> Dict[str, float]:
    out: Dict[str, float] = {}
    for key, pat in PATTERNS.items():
        m = pat.search(text)
        if not m:
            raise RuntimeError(f"Could not find '{key}' in output:\n{text}")
        val_str = m.group(1).replace(",", "")
        out[key] = float(val_str)
    return out

# ╭──────────────────────────────────────────────────────────╮
# │ Run the benchmark once, stream output, capture for parse │
# ╰──────────────────────────────────────────────────────────╯
def run_once(script: str) -> Dict[str, float]:
    env = os.environ.copy()
    env["PYTHONUNBUFFERED"] = "1"        

    proc = subprocess.Popen(
        ["python", script],
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,
        bufsize=0,                     
        env=env,
    )

    captured = io.StringIO()
    assert proc.stdout is not None

    with proc.stdout:
        while True:
            chunk = proc.stdout.read(4096)
            if not chunk:
                break
            text = chunk.decode(errors="replace")
            sys.stdout.write(text)
            sys.stdout.flush()
            captured.write(text)

    proc.wait()
    if proc.returncode != 0:
        raise RuntimeError(f"{script} exited with code {proc.returncode}")

    return parse_stdout(captured.getvalue())

# ╭──────────────────────────────────────────────────────────╮
# │ Aggregate metric lists with arithmetic means             │
# ╰──────────────────────────────────────────────────────────╯
def aggregate(results: List[Dict[str, float]]) -> Dict[str, float]:
    keys = results[0].keys()
    return {k: stats.mean(r[k] for r in results) for k in keys}

# ╭──────────────────────────────────────────────────────────╮
# │ Driver                                                   │
# ╰──────────────────────────────────────────────────────────╯
def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("-n", "--runs", type=int, default=3,
                        help="Number of times to run the benchmark (default: 3)")
    parser.add_argument("-s", "--script", default="hsr_speed_test.py",
                        help="Path to benchmark script")
    parser.add_argument("-o", "--output",
                        help="CSV file for per‑run stats and averages")
    args = parser.parse_args()

    all_results: List[Dict[str, float]] = []
    for idx in range(1, args.runs + 1):
        print(f"\n=== Run {idx}/{args.runs} ===")
        all_results.append(run_once(args.script))

    ave = aggregate(all_results)

    # ── Pretty‑print averages (files metric omitted) ──────────────
    LABELS = {
        "molecules":    "Molecules processed         ",
        "avg_fp":       "Average fp time (s)         ",
        "max_fp":       "Max fp time (s)             ",
        "min_fp":       "Min fp time (s)             ",
        "seq_fp_time":  "Sequential fp time (s)      ",
        "sim_count":    "Similarity comps (#)        ",
        "sim_time":     "Similarity comps time (s)   ",
        "wall_fp":      "Wall‑clock fp phase (s)     ",
        "wall_dl":      "Wall‑clock dl phase (s)     ",
        "runtime":      "End‑to‑end runtime (s)      ",
    }

    print("\n────────── AVERAGED RESULTS ──────────")
    for key in LABELS:
        val = ave[key]
        # print similarity count as int, others as float
        if key == "sim_count" or key == "molecules":
            print(f"{LABELS[key]}: {val:,.0f}")
        else:
            print(f"{LABELS[key]}: {val:,.4f}")

    # ── Optional CSV export ───────────────────────────────────────
    if args.output:
        path = Path(args.output)
        with path.open("w", newline="") as fh:
            writer = csv.writer(fh)
            header = ["run"] + list(all_results[0].keys())
            writer.writerow(header)
            for idx, row in enumerate(all_results, 1):
                writer.writerow([idx] + [row[k] for k in header[1:]])
            writer.writerow(["average"] + [ave[k] for k in header[1:]])
        print(f"\n📄 Per‑run stats + averages saved to {path.resolve()}")

if __name__ == "__main__":
    main()
