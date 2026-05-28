#!/usr/bin/env python3
"""Parse STAP++ Q4 verification outputs and write CSV summary tables."""
from __future__ import annotations

import argparse
import csv
import json
import re
from pathlib import Path
from statistics import mean
from typing import Dict, List, Tuple

FLOAT = r"[-+]?\d*\.\d+e[-+]?\d+|[-+]?\d+\.\d*e[-+]?\d+|[-+]?\d+e[-+]?\d+"


def parse_output(path: Path) -> Tuple[Dict[int, Tuple[float, float, float]], Dict[int, Tuple[float, float, float]]]:
    displacements: Dict[int, Tuple[float, float, float]] = {}
    stresses: Dict[int, Tuple[float, float, float]] = {}
    section = None
    for raw in path.read_text(errors="ignore").splitlines():
        line = raw.rstrip()
        if "D I S P L A C E M E N T S" in line:
            section = "disp"
            continue
        if "S T R E S S" in line:
            section = "stress"
            continue
        if "S O L U T I O N" in line:
            section = None
        m = re.match(rf"^\s*(\d+)\s+({FLOAT})\s+({FLOAT})\s+({FLOAT})\s*$", line, re.IGNORECASE)
        if not m:
            continue
        idx = int(m.group(1))
        vals = (float(m.group(2)), float(m.group(3)), float(m.group(4)))
        if section == "disp":
            displacements[idx] = vals
        elif section == "stress":
            stresses[idx] = vals
    return displacements, stresses


def write_patch_summary(case: Dict, out_path: Path, rows_path: Path) -> None:
    disp, stress = parse_output(out_path)
    expected_nodes = case["expected_nodes"]
    node_rows: List[Dict] = []
    max_ux = 0.0
    max_uy = 0.0
    for nid_s, exp in expected_nodes.items():
        nid = int(nid_s)
        got = disp[nid]
        err_ux = got[0] - exp["ux"]
        err_uy = got[1] - exp["uy"]
        max_ux = max(max_ux, abs(err_ux))
        max_uy = max(max_uy, abs(err_uy))
        node_rows.append({
            "node": nid,
            "ux": got[0],
            "ux_expected": exp["ux"],
            "ux_error": err_ux,
            "uy": got[1],
            "uy_expected": exp["uy"],
            "uy_error": err_uy,
        })
    rows_path.parent.mkdir(parents=True, exist_ok=True)
    with rows_path.open("w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=list(node_rows[0].keys()))
        writer.writeheader()
        writer.writerows(node_rows)

    expected_stress = case["expected_stress"]
    sx = mean(v[0] for v in stress.values())
    sy = mean(v[1] for v in stress.values())
    txy = mean(v[2] for v in stress.values())
    summary = [{
        "case": case["name"],
        "max_abs_ux_error": max_ux,
        "max_abs_uy_error": max_uy,
        "sigma_x_avg": sx,
        "sigma_x_expected": expected_stress["sigma_x"],
        "sigma_x_error": sx - expected_stress["sigma_x"],
        "sigma_y_avg": sy,
        "sigma_y_expected": expected_stress["sigma_y"],
        "sigma_y_error": sy - expected_stress["sigma_y"],
        "tau_xy_avg": txy,
        "tau_xy_expected": expected_stress["tau_xy"],
        "tau_xy_error": txy - expected_stress["tau_xy"],
    }]
    with (rows_path.parent / "q4_patch_summary.csv").open("w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=list(summary[0].keys()))
        writer.writeheader()
        writer.writerows(summary)


def write_convergence_summary(cases: List[Dict], dat_dir: Path, table_dir: Path) -> None:
    rows: List[Dict] = []
    for case in cases:
        out_path = dat_dir / f"{case['name']}.out"
        disp, stress = parse_output(out_path)
        # Right edge has maximum node displacement for these left-fixed cases.
        max_ux = max(v[0] for v in disp.values())
        max_uy = max(v[1] for v in disp.values())
        min_uy = min(v[1] for v in disp.values())
        rows.append({
            "case": case["name"],
            "family": case["family"],
            "n": case["n"],
            "max_ux": max_ux,
            "max_uy": max_uy,
            "min_uy": min_uy,
            "avg_sigma_x": mean(v[0] for v in stress.values()),
            "avg_sigma_y": mean(v[1] for v in stress.values()),
            "avg_tau_xy": mean(v[2] for v in stress.values()),
        })
    # Add self-convergence deltas against the finest mesh in each family.
    finest = {fam: max([r for r in rows if r["family"] == fam], key=lambda r: r["n"]) for fam in {r["family"] for r in rows}}
    for r in rows:
        ref = finest[r["family"]]
        metric = "max_ux" if r["family"] == "tension" else "max_uy"
        ref_val = ref[metric]
        r["reference_metric"] = metric
        r["reference_value_8x8"] = ref_val
        r["abs_diff_to_8x8"] = abs(r[metric] - ref_val)
        r["rel_diff_to_8x8"] = abs(r[metric] - ref_val) / abs(ref_val) if abs(ref_val) > 0 else 0.0
    table_dir.mkdir(parents=True, exist_ok=True)
    with (table_dir / "q4_convergence_summary.csv").open("w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=list(rows[0].keys()))
        writer.writeheader()
        writer.writerows(sorted(rows, key=lambda r: (r["family"], r["n"])))


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--manifest", default="data/q4_generated/manifest.json")
    parser.add_argument("--dat-dir", default="data/q4_generated")
    parser.add_argument("--table-dir", default="reports/tables")
    args = parser.parse_args()
    manifest = json.loads(Path(args.manifest).read_text(encoding="utf-8"))
    dat_dir = Path(args.dat_dir)
    table_dir = Path(args.table_dir)
    patch = next(c for c in manifest["cases"] if c["type"] == "patch")
    write_patch_summary(patch, dat_dir / f"{patch['name']}.out", table_dir / "q4_patch_nodes.csv")
    conv = [c for c in manifest["cases"] if c["type"] == "convergence"]
    write_convergence_summary(conv, dat_dir, table_dir)
    print(f"Wrote CSV tables to {table_dir}")


if __name__ == "__main__":
    main()
