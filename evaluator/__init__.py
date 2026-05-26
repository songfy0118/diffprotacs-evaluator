"""
DiffPROTACs Evaluator
=====================
A drop-in evaluation pipeline for DiffPROTACs: turns 3D `.xyz` outputs
into SMILES, scores them on validity / QED / SA / Lipinski / alerts /
novelty, and optionally runs DockStream docking.

Quick start:
    from evaluator.scorers import default_scorer
    from evaluator.xyz_to_smiles import xyz_to_smiles

    smi = xyz_to_smiles("path/to/some.xyz")
    scorer = default_scorer()
    print(scorer.score_one(smi))

Or run end-to-end:
    python -m evaluator.evaluate --xyz_dir <DiffPROTACs output> --out eval.csv
"""

from .xyz_to_smiles import xyz_to_mol, mol_to_smiles, xyz_to_smiles, collect_xyz_files
from .scorers import (
    CompositeScorer,
    ValidityScorer,
    QEDScorer,
    SAScorer,
    LipinskiScorer,
    LinkerLengthScorer,
    SubstructureAlerts,
    NoveltyScorer,
    default_scorer,
)

__all__ = [
    "xyz_to_mol",
    "mol_to_smiles",
    "xyz_to_smiles",
    "collect_xyz_files",
    "CompositeScorer",
    "ValidityScorer",
    "QEDScorer",
    "SAScorer",
    "LipinskiScorer",
    "LinkerLengthScorer",
    "SubstructureAlerts",
    "NoveltyScorer",
    "default_scorer",
]
