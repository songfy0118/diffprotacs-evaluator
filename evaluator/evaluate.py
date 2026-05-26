"""
evaluate.py
===========
End-to-end evaluator for DiffPROTACs outputs.

Usage
-----
    python -m evaluator.evaluate \
        --xyz_dir   path/to/diffprotacs/output_dir   \
        --out       results.csv                      \
        --reference reference_smiles.txt   (optional) \
        --docking_config configs/dockstream_template.json  (optional) \
        --dockstream_root path/to/DockStream-master-master  (optional)

What it does
------------
1. Walk the DiffPROTACs output directory and collect sampled .xyz files
   (skipping `true`/`frag` reference files).
2. Convert each xyz -> RDKit Mol -> canonical SMILES.
3. Score every SMILES with the default RDKit-based CompositeScorer
   (validity, QED, SA, Lipinski, size window, substructure alerts,
   optional novelty vs reference).
4. If a DockStream config is provided, batch all VALID SMILES through
   DockStream/Vina and append `raw_dock` + `dock_score` columns.
5. Write per-molecule results to CSV and print a summary block
   (validity rate, uniqueness, mean QED, mean docking, top-K table).

Output CSV columns
------------------
    uuid, sample_idx, xyz_path, smiles, validity, qed, sa_proxy,
    lipinski_fraction, size_window, alerts_pass, novelty, total_score,
    raw_dock, dock_score   (last two only if docking enabled)

Exit code 0 on success. Errors during scoring of individual molecules
are caught and logged; only catastrophic failures (e.g. bad paths,
malformed CSV) raise.
"""

from __future__ import annotations

import argparse
import csv
import os
import sys
from collections import Counter
from pathlib import Path
from typing import Dict, List, Optional

# Allow running both as module (-m evaluator.evaluate) and as script
try:
    from .xyz_to_smiles import collect_xyz_files, xyz_to_mol, mol_to_smiles
    from .scorers import default_scorer, NoveltyScorer
    from .docking_bridge import DockStreamScorer
except ImportError:
    sys.path.insert(0, str(Path(__file__).parent.parent))
    from evaluator.xyz_to_smiles import collect_xyz_files, xyz_to_mol, mol_to_smiles
    from evaluator.scorers import default_scorer, NoveltyScorer
    from evaluator.docking_bridge import DockStreamScorer


def _parse_uuid_idx(xyz_path: str) -> Dict[str, str]:
    """
    DiffPROTACs writes `<output_dir>/<uuid>/<sample_idx>_.xyz`.
    Extract uuid and sample_idx for the report.
    """
    p = Path(xyz_path)
    sample_idx = p.stem  # e.g. "0_" or "0"
    uuid = p.parent.name
    return {"uuid": uuid, "sample_idx": sample_idx}


def _load_reference(path: Optional[str]) -> List[str]:
    if path is None:
        return []
    if not os.path.exists(path):
        print(f"[evaluate] reference file not found: {path} — skipping novelty")
        return []
    smis = []
    with open(path) as f:
        for line in f:
            s = line.strip().split()[0] if line.strip() else ""
            if s:
                smis.append(s)
    print(f"[evaluate] loaded {len(smis)} reference SMILES")
    return smis


def run(args: argparse.Namespace) -> int:
    # 1. Collect files
    xyz_paths = collect_xyz_files(args.xyz_dir)
    if not xyz_paths:
        print(f"[evaluate] no .xyz files found under {args.xyz_dir}")
        return 1
    print(f"[evaluate] found {len(xyz_paths)} sampled molecules")

    # 2. xyz -> SMILES
    print("[evaluate] converting xyz -> SMILES …")
    rows: List[Dict] = []
    n_valid = 0
    for i, xp in enumerate(xyz_paths):
        meta = _parse_uuid_idx(xp)
        mol = xyz_to_mol(xp)
        smi = mol_to_smiles(mol)
        if smi:
            n_valid += 1
        rows.append({**meta, "xyz_path": xp, "smiles": smi or ""})
        if (i + 1) % 100 == 0:
            print(f"   {i + 1}/{len(xyz_paths)} processed ({n_valid} valid so far)")

    print(f"[evaluate] validity rate (xyz -> mol): {n_valid}/{len(xyz_paths)} "
          f"= {n_valid / len(xyz_paths):.2%}")

    # 3. RDKit scoring
    print("[evaluate] computing RDKit scores …")
    ref_smiles = _load_reference(args.reference)
    scorer = default_scorer(reference_smiles=ref_smiles)
    for r in rows:
        scores = scorer.score_one(r["smiles"] or None)
        r.update(scores)

    # 4. (optional) Docking
    if args.docking_config and args.dockstream_root:
        print("[evaluate] running DockStream docking on valid SMILES …")
        valid_smis = [r["smiles"] for r in rows if r["smiles"]]
        # de-duplicate to avoid redocking the same molecule
        unique_smis = sorted(set(valid_smis))
        print(f"   {len(unique_smis)} unique SMILES to dock")
        try:
            ds = DockStreamScorer(
                dockstream_root=args.dockstream_root,
                config_path=args.docking_config,
                python_executable=args.dockstream_python,
                timeout_sec=args.dockstream_timeout,
            )
            dock_results = ds.dock(unique_smis)
        except FileNotFoundError as e:
            print(f"[evaluate] docking skipped: {e}")
            dock_results = {}
        for r in rows:
            d = dock_results.get(r["smiles"], {"raw_dock": "", "dock_score": ""})
            r["raw_dock"] = d["raw_dock"]
            r["dock_score"] = d["dock_score"]
    elif args.docking_config or args.dockstream_root:
        print("[evaluate] need BOTH --docking_config and --dockstream_root for docking; skipping")

    # 5. Write CSV
    print(f"[evaluate] writing results to {args.out}")
    # Determine union of keys, but force a sensible column order
    preferred = [
        "uuid", "sample_idx", "xyz_path", "smiles",
        "validity", "qed", "sa_proxy", "lipinski_fraction",
        "size_window", "alerts_pass", "novelty", "total_score",
        "raw_dock", "dock_score",
    ]
    all_keys = set()
    for r in rows:
        all_keys.update(r.keys())
    fieldnames = [k for k in preferred if k in all_keys] + sorted(all_keys - set(preferred))

    with open(args.out, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=fieldnames)
        w.writeheader()
        for r in rows:
            w.writerow(r)

    # 6. Summary
    _print_summary(rows)
    return 0


def _print_summary(rows: List[Dict]) -> None:
    n = len(rows)
    valid = [r for r in rows if r["smiles"]]
    n_valid = len(valid)
    n_unique = len({r["smiles"] for r in valid})
    print()
    print("=" * 60)
    print("DiffPROTACs Evaluation Summary")
    print("=" * 60)
    print(f"  Total samples         : {n}")
    print(f"  Valid (parses as mol) : {n_valid}  ({n_valid / max(1, n):.2%})")
    print(f"  Unique SMILES         : {n_unique}  "
          f"({n_unique / max(1, n_valid):.2%} of valid)")

    def _mean(key):
        vs = [float(r[key]) for r in valid if r.get(key) not in ("", None)]
        return sum(vs) / len(vs) if vs else float("nan")

    for key in ["qed", "sa_proxy", "lipinski_fraction", "novelty",
                "alerts_pass", "total_score", "raw_dock", "dock_score"]:
        if any(r.get(key) not in ("", None) for r in valid):
            print(f"  Mean {key:<20}: {_mean(key):.4f}")

    # Top-5 by total_score
    if valid:
        top = sorted(valid, key=lambda r: -float(r.get("total_score") or 0))[:5]
        print()
        print("  Top-5 by total_score:")
        for i, r in enumerate(top, 1):
            print(f"    {i}. score={float(r['total_score']):.3f}  smi={r['smiles'][:80]}")
    print("=" * 60)


def main():
    parser = argparse.ArgumentParser(
        description="Evaluate molecules generated by DiffPROTACs.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument(
        "--xyz_dir", required=True,
        help="DiffPROTACs sample.py output directory "
             "(contains <uuid>/<i>_.xyz files).",
    )
    parser.add_argument(
        "--out", default="diffprotacs_eval.csv",
        help="Where to write the per-molecule evaluation CSV.",
    )
    parser.add_argument(
        "--reference", default=None,
        help="Optional path to a .smi/.txt file of reference SMILES (one per "
             "line). Used for novelty scoring (e.g. the training set, or "
             "PROTAC-DB). Without this, novelty is 1.0 for every valid mol.",
    )
    # Docking opts
    parser.add_argument(
        "--docking_config", default=None,
        help="Path to a DockStream JSON config (e.g. configs/dockstream_"
             "template.json filled in for your target). Enables docking.",
    )
    parser.add_argument(
        "--dockstream_root", default=None,
        help="Path to your DockStream-master-master folder. Required if "
             "--docking_config is set.",
    )
    parser.add_argument(
        "--dockstream_python", default="python",
        help="Interpreter that has DockStream's dependencies installed.",
    )
    parser.add_argument(
        "--dockstream_timeout", type=int, default=1800,
        help="Hard timeout (seconds) for the docking subprocess.",
    )
    args = parser.parse_args()
    sys.exit(run(args))


if __name__ == "__main__":
    main()
