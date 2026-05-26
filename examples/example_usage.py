"""
example_usage.py
================
Three quick examples of how to use the DiffPROTACs evaluator from Python
(rather than the CLI). Run this from the repo root:

    python examples/example_usage.py
"""

import os
import sys

# Make the repo root importable when this script is run directly
HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
if ROOT not in sys.path:
    sys.path.insert(0, ROOT)

# Also chdir to repo root so the demo_data/ relative paths resolve
os.chdir(ROOT)

from evaluator import (
    xyz_to_smiles,
    xyz_to_mol,
    default_scorer,
    CompositeScorer,
    QEDScorer,
    SAScorer,
    SubstructureAlerts,
    NoveltyScorer,
    collect_xyz_files,
)


# -----------------------------------------------------------------------------
# Example 1: convert a single xyz to SMILES
# -----------------------------------------------------------------------------
def example_1_single_file():
    print("\n[Example 1] Single .xyz -> SMILES")
    print("-" * 40)
    smi = xyz_to_smiles("demo_data/demo_output/amide_linker1/0_.xyz")
    print(f"SMILES: {smi}")


# -----------------------------------------------------------------------------
# Example 2: score a SMILES with the default scorer
# -----------------------------------------------------------------------------
def example_2_default_scorer():
    print("\n[Example 2] Default scorer")
    print("-" * 40)
    scorer = default_scorer()
    result = scorer.score_one("CC(=O)NCCCCCNC(C)=O")  # the amide linker
    for k, v in result.items():
        print(f"  {k:<20} : {v:.4f}")


# -----------------------------------------------------------------------------
# Example 3: build a CUSTOM scorer combination
# -----------------------------------------------------------------------------
def example_3_custom_scorer():
    """
    Suppose you only care about QED + tox alerts + novelty vs your
    PROTAC-DB 3.0 reference set, with extra weight on alerts (hard gate).
    """
    print("\n[Example 3] Custom scorer composition")
    print("-" * 40)
    reference_smiles = [
        # In practice load these from PROTAC-DB; here is a tiny sample.
        "CC(=O)NCCCCCNC(C)=O",
        "COCCOCCOCCOC",
    ]
    scorer = CompositeScorer(
        components=[
            (QEDScorer(), 1.0),
            (SAScorer(), 0.5),
            (SubstructureAlerts(), 2.0),  # double-weighted hard gate
            (NoveltyScorer(reference_smiles), 1.0),
        ],
        mode="product",
    )
    for smi in [
        "c1ccc2ccccc2c1",              # naphthalene
        "CC(=O)NCCCCCNC(C)=O",         # IN the reference -> novelty=0 -> total=0
        "Cc1ccc(NC(=O)CCNC(=O)C)cc1",  # novel amide
    ]:
        r = scorer.score_one(smi)
        print(f"  {smi[:40]:<40}  total={r['total_score']:.3f}  "
              f"novelty={r['novelty']:.2f}  qed={r['qed']:.2f}")


# -----------------------------------------------------------------------------
# Example 4: batch over a DiffPROTACs output directory
# -----------------------------------------------------------------------------
def example_4_batch():
    print("\n[Example 4] Batch over a sample directory")
    print("-" * 40)
    paths = collect_xyz_files("demo_data/demo_output")
    print(f"  found {len(paths)} sample files")
    scorer = default_scorer()
    n_valid = 0
    qeds = []
    for p in paths[:8]:
        smi = xyz_to_smiles(p)
        if smi:
            n_valid += 1
            qeds.append(scorer.score_one(smi)["qed"])
    print(f"  valid: {n_valid}/{min(8, len(paths))}")
    if qeds:
        print(f"  mean QED: {sum(qeds) / len(qeds):.3f}")


if __name__ == "__main__":
    example_1_single_file()
    example_2_default_scorer()
    example_3_custom_scorer()
    example_4_batch()
