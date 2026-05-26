"""
docking_bridge.py
=================
Optional bridge to DockStream for AutoDock Vina (or other backends)
docking scores.

This module is **optional**. The core evaluator runs without it. We use
DockStream the way Protac-invent does — via the JSON config + CLI route,
so we don't import its python modules and avoid env conflicts.

Workflow:
    smiles list -> write smiles file -> call DockStream's docker.py
                -> parse the produced scores CSV -> return dict[smiles]=score

Vina scores are NEGATIVE (kcal/mol). Lower (more negative) = better
binding. We expose two values:
    raw_dock   : original Vina score, e.g. -8.4
    dock_score : mapped to [0,1] via sigmoid, higher = better
                 (default: midpoint at -7, slope 1.0)

If DockStream isn't installed / configured, calls return empty dicts
and a warning. The CompositeScorer handles missing docking gracefully.

Typical config file: configs/dockstream_template.json — fill in the
binary_location, receptor_pdbqt_path, and search_space for your target.
"""

from __future__ import annotations

import os
import csv
import json
import math
import shutil
import subprocess
import tempfile
from pathlib import Path
from typing import Dict, List, Optional


def _sigmoid(x: float, k: float = 1.0, x0: float = 0.0) -> float:
    try:
        return 1.0 / (1.0 + math.exp(-k * (x - x0)))
    except OverflowError:
        return 0.0 if (x - x0) < 0 else 1.0


def vina_to_unit(raw: float, midpoint: float = -7.0, slope: float = 1.0) -> float:
    """
    Map Vina score (more negative = better, typically -6 to -12) to [0,1]
    where higher = better. Default midpoint -7 kcal/mol corresponds to
    "OK" binding -> ~0.5. -10 -> ~0.95. -4 -> ~0.05.
    """
    # We invert sign so larger sigmoid arg means better
    return _sigmoid(-raw, k=slope, x0=-midpoint)


class DockStreamScorer:
    """
    Thin wrapper around DockStream's `docker.py` CLI.

    Parameters
    ----------
    dockstream_root : path to the DockStream-master-master folder
    config_path     : a docking config JSON (see DockStream's docs and the
                      template under configs/dockstream_template.json)
    python_executable : interpreter that has DockStream deps installed.
                        Default 'python' assumes you're already inside
                        DockStream's conda env.
    """

    name = "docking"

    def __init__(
        self,
        dockstream_root: str,
        config_path: str,
        python_executable: str = "python",
        timeout_sec: int = 1800,
        midpoint: float = -7.0,
        slope: float = 1.0,
    ):
        self.root = Path(dockstream_root)
        self.config_path = Path(config_path)
        self.python = python_executable
        self.timeout = timeout_sec
        self.midpoint = midpoint
        self.slope = slope

        if not self.root.exists():
            raise FileNotFoundError(f"DockStream root not found: {self.root}")
        if not self.config_path.exists():
            raise FileNotFoundError(f"Docking config not found: {self.config_path}")
        if not (self.root / "docker.py").exists():
            raise FileNotFoundError(
                f"docker.py not found under {self.root}; check dockstream_root"
            )

    def dock(self, smiles: List[str]) -> Dict[str, Dict[str, float]]:
        """
        Dock a list of SMILES. Returns {smiles: {'raw_dock': float,
        'dock_score': float}}. SMILES that fail to dock get a 'raw_dock'
        of +inf and 'dock_score' of 0.0.
        """
        if not smiles:
            return {}

        # Each call gets its own temp output dir so concurrent runs don't
        # stomp on each other. We rewrite the config file's score path to
        # point into this dir.
        out_dir = Path(tempfile.mkdtemp(prefix="dockstream_"))
        scores_csv = out_dir / "scores.csv"
        poses_sdf = out_dir / "poses.sdf"
        patched_cfg = out_dir / "config.json"
        try:
            with open(self.config_path) as f:
                cfg = json.load(f)
            # Best-effort patch: redirect output paths into our tmp dir.
            # The exact JSON keys depend on the DockStream config schema;
            # this matches the configs shipped under Protac-invent.
            try:
                runs = cfg["docking"]["docking_runs"]
                for run in runs:
                    run["output"]["scores"]["scores_path"] = str(scores_csv)
                    run["output"]["poses"]["poses_path"] = str(poses_sdf)
            except KeyError:
                pass  # leave config unchanged; user controls output
            with open(patched_cfg, "w") as f:
                json.dump(cfg, f, indent=2)

            cmd = [
                self.python,
                str(self.root / "docker.py"),
                "-conf",
                str(patched_cfg),
                "-smiles",
                ";".join(smiles),
                "-print_scores",
            ]
            try:
                proc = subprocess.run(
                    cmd, capture_output=True, text=True, timeout=self.timeout
                )
            except subprocess.TimeoutExpired:
                print(f"[DockStreamScorer] timed out after {self.timeout}s")
                return {s: {"raw_dock": float("inf"), "dock_score": 0.0} for s in smiles}

            if proc.returncode != 0:
                print(f"[DockStreamScorer] DockStream failed (rc={proc.returncode}):")
                print(proc.stderr[-2000:])
                return {s: {"raw_dock": float("inf"), "dock_score": 0.0} for s in smiles}

            return self._parse_scores(smiles, scores_csv, proc.stdout)
        finally:
            shutil.rmtree(out_dir, ignore_errors=True)

    def _parse_scores(
        self, smiles: List[str], scores_csv: Path, stdout: str
    ) -> Dict[str, Dict[str, float]]:
        """Read scores from the CSV DockStream writes; fall back to stdout."""
        raw_by_smi: Dict[str, float] = {s: float("inf") for s in smiles}

        # Preferred path: CSV
        if scores_csv.exists():
            try:
                with open(scores_csv) as f:
                    reader = csv.DictReader(f)
                    for row in reader:
                        smi = row.get("smiles") or row.get("SMILES")
                        score_str = (
                            row.get("score") or row.get("Score") or row.get("ADV") or ""
                        )
                        if smi and score_str:
                            try:
                                raw_by_smi[smi] = float(score_str)
                            except ValueError:
                                pass
            except Exception as e:
                print(f"[DockStreamScorer] CSV parse failed: {e}")

        # Fallback: try to read print_scores stdout
        if all(v == float("inf") for v in raw_by_smi.values()):
            for line in stdout.splitlines():
                # DockStream with -print_scores writes: "SMILES\tscore"
                if "\t" in line:
                    parts = line.split("\t")
                    if len(parts) >= 2:
                        try:
                            raw_by_smi[parts[0]] = float(parts[1])
                        except ValueError:
                            continue

        return {
            s: {
                "raw_dock": raw_by_smi[s],
                "dock_score": vina_to_unit(
                    raw_by_smi[s], midpoint=self.midpoint, slope=self.slope
                )
                if raw_by_smi[s] != float("inf")
                else 0.0,
            }
            for s in smiles
        }

    # Adapter for CompositeScorer (per-molecule API).
    # In practice you batch all SMILES through `.dock(list)` for speed;
    # this method exists so docking can also be dropped into a Composite.
    def score(self, smiles: Optional[str]) -> float:
        if not smiles:
            return 0.0
        result = self.dock([smiles]).get(smiles)
        return result["dock_score"] if result else 0.0
