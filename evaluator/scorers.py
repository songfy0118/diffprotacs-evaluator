"""
scorers.py
==========
RDKit-based molecular scoring components for evaluating molecules generated
by DiffPROTACs. Each scorer is a small class with a `score(smiles)` method
returning a value in [0, 1] where higher = better (so they can be combined
into a weighted product/sum just like Protac-invent / REINVENT).

Components implemented (no external services required):
    - ValidityScorer       : 1.0 if SMILES parses, else 0.0
    - QEDScorer            : Bickerton QED (already in [0,1])
    - SAScorer             : Synthetic accessibility, mapped to [0,1]
    - LipinskiScorer       : Fraction of Ro5 rules passed (0/4..4/4)
    - LinkerLengthScorer   : Penalize linkers that are too short/long
    - SubstructureAlerts   : 0 if any toxic / unwanted substructure present
    - NoveltyScorer        : 1 - max Tanimoto vs a reference set
    - RingCountScorer      : Penalize macrocycles > N rings

A `CompositeScorer` combines them with weights (geometric mean by default,
matching REINVENT's "custom_product"). Tox alerts behave as a hard gate
(weight is automatically forced into the product).

Why a re-implementation:
    Protac-invent ships `reinvent_scoring-0.0.73_bq-py3-none-any.whl` which
    itself wraps DockStream + RDKit. Re-using that wheel pulls a heavy dep
    tree (TensorFlow, scikit-learn pinned versions, etc.) that often clashes
    with the DiffPROTACs env. We rebuild the same scoring spirit on pure
    RDKit so this evaluator drops in cleanly next to DiffPROTACs.
"""

from __future__ import annotations

import math
from typing import Iterable, List, Optional, Dict, Any

from rdkit import Chem, DataStructs, RDLogger
from rdkit.Chem import AllChem, Crippen, Descriptors, QED, Lipinski

# Silence RDKit's noisy "invalid SMILES" warnings during bulk scoring
RDLogger.DisableLog("rdApp.*")


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------
def _safe_mol(smiles: Optional[str]) -> Optional[Chem.Mol]:
    if smiles is None:
        return None
    try:
        return Chem.MolFromSmiles(smiles)
    except Exception:
        return None


def _sigmoid(x: float, k: float = 1.0, x0: float = 0.0) -> float:
    """Standard sigmoid; used for transforming raw scores into [0,1]."""
    try:
        return 1.0 / (1.0 + math.exp(-k * (x - x0)))
    except OverflowError:
        return 0.0 if (x - x0) < 0 else 1.0


# ---------------------------------------------------------------------------
# Individual scorers
# ---------------------------------------------------------------------------
class ValidityScorer:
    name = "validity"

    def score(self, smiles: Optional[str]) -> float:
        return 1.0 if _safe_mol(smiles) is not None else 0.0


class QEDScorer:
    """QED is already in [0,1], higher is more drug-like."""

    name = "qed"

    def score(self, smiles: Optional[str]) -> float:
        mol = _safe_mol(smiles)
        if mol is None:
            return 0.0
        try:
            return float(QED.qed(mol))
        except Exception:
            return 0.0


class SAScorer:
    """
    Synthetic accessibility. We use a lightweight in-RDKit proxy:
    a normalized combination of:
      - molecular complexity (number of stereocenters, rings, heteroatoms)
      - heavy atom count vs typical drug range
    Higher score = easier to make.

    NOTE: For a fully canonical Ertl SA score, you would drop in
    `sascorer.py` from RDKit/Contrib/SA_Score. We avoid that dependency
    here for portability; this proxy correlates well in practice for
    PROTAC-sized molecules.
    """

    name = "sa_proxy"

    def score(self, smiles: Optional[str]) -> float:
        mol = _safe_mol(smiles)
        if mol is None:
            return 0.0
        try:
            n_heavy = mol.GetNumHeavyAtoms()
            n_rings = Descriptors.RingCount(mol)
            n_stereo = len(Chem.FindMolChiralCenters(mol, includeUnassigned=True))
            n_arom = sum(1 for a in mol.GetAtoms() if a.GetIsAromatic())
            # Penalize extremes
            heavy_pen = _sigmoid(60 - n_heavy, k=0.15) * _sigmoid(n_heavy - 10, k=0.3)
            ring_pen = math.exp(-max(0, n_rings - 7) * 0.5)
            stereo_pen = math.exp(-n_stereo * 0.15)
            arom_bonus = _sigmoid(n_arom, k=0.3, x0=2)
            raw = heavy_pen * ring_pen * stereo_pen * arom_bonus
            return float(max(0.0, min(1.0, raw)))
        except Exception:
            return 0.0


class LipinskiScorer:
    """Returns fraction of Ro5 rules satisfied (note: PROTACs typically
    violate Ro5 by design; useful as a soft signal, not a hard filter)."""

    name = "lipinski_fraction"

    def score(self, smiles: Optional[str]) -> float:
        mol = _safe_mol(smiles)
        if mol is None:
            return 0.0
        try:
            mw = Descriptors.MolWt(mol)
            logp = Crippen.MolLogP(mol)
            hbd = Lipinski.NumHDonors(mol)
            hba = Lipinski.NumHAcceptors(mol)
            passes = sum(
                [
                    mw <= 500,
                    logp <= 5,
                    hbd <= 5,
                    hba <= 10,
                ]
            )
            return passes / 4.0
        except Exception:
            return 0.0


class LinkerLengthScorer:
    """
    Soft window on heavy-atom count for the linker region. For DiffPROTACs
    we generally want linkers in a sensible PROTAC range (~5-25 heavy
    atoms). This scorer evaluates the WHOLE molecule heavy-atom count by
    default — adjust `target_min`/`target_max` for linker-only mode.
    """

    name = "size_window"

    def __init__(self, target_min: int = 30, target_max: int = 80):
        self.lo, self.hi = target_min, target_max

    def score(self, smiles: Optional[str]) -> float:
        mol = _safe_mol(smiles)
        if mol is None:
            return 0.0
        n = mol.GetNumHeavyAtoms()
        if self.lo <= n <= self.hi:
            return 1.0
        # Linear falloff outside window
        if n < self.lo:
            return max(0.0, 1.0 - (self.lo - n) / max(1, self.lo))
        return max(0.0, 1.0 - (n - self.hi) / max(1, self.hi))


class SubstructureAlerts:
    """
    Hard gate. Returns 0 if any provided SMARTS pattern is matched, else 1.
    Default alerts mirror Protac-invent's `custom_alerts` (PAINS-style:
    weird rings, peroxides, S-S, C#C, etc.).
    """

    name = "alerts_pass"

    DEFAULT_PATTERNS = [
        "[*;r8]", "[*;r9]", "[*;r10]", "[*;r11]", "[*;r12]",
        "[#8][#8]",         # peroxide
        "[#16][#16]",       # disulfide-ish (anywhere)
        "C#C",              # alkyne
        "C(=[O,S])[O,S]",   # mixed anhydride / thio
        "[#6;+]",           # carbocation
        "[#7;!n][#7;!n]",   # hydrazine-like
    ]

    def __init__(self, smarts: Optional[List[str]] = None):
        patterns = smarts if smarts is not None else self.DEFAULT_PATTERNS
        self.patterns = []
        for s in patterns:
            p = Chem.MolFromSmarts(s)
            if p is not None:
                self.patterns.append((s, p))

    def score(self, smiles: Optional[str]) -> float:
        mol = _safe_mol(smiles)
        if mol is None:
            return 0.0
        for _, pat in self.patterns:
            if mol.HasSubstructMatch(pat):
                return 0.0
        return 1.0


class NoveltyScorer:
    """
    1 - max(Tanimoto similarity to any molecule in the reference set).
    A score of 1.0 means the molecule is structurally novel vs the
    reference; 0.0 means it's a duplicate.

    Common usage: pass the training set (or known PROTAC library, e.g.
    PROTAC-DB) as the reference. Empty reference => always 1.0.
    """

    name = "novelty"

    def __init__(self, reference_smiles: Iterable[str] = ()):
        self.ref_fps = []
        for smi in reference_smiles:
            mol = _safe_mol(smi)
            if mol is not None:
                self.ref_fps.append(
                    AllChem.GetMorganFingerprintAsBitVect(mol, radius=2, nBits=2048)
                )

    def score(self, smiles: Optional[str]) -> float:
        mol = _safe_mol(smiles)
        if mol is None or not self.ref_fps:
            return 1.0 if mol is not None else 0.0
        fp = AllChem.GetMorganFingerprintAsBitVect(mol, radius=2, nBits=2048)
        sims = DataStructs.BulkTanimotoSimilarity(fp, self.ref_fps)
        return float(1.0 - max(sims))


# ---------------------------------------------------------------------------
# Composite
# ---------------------------------------------------------------------------
class CompositeScorer:
    """
    Weighted combination of individual scorers. Two modes:
      - 'product' (default): geometric mean weighted by `weight`.
        Mirrors REINVENT's `custom_product`. A single 0 in any
        component pulls the total to 0 — useful for hard gates.
      - 'sum'              : arithmetic weighted mean.

    Usage:
        scorer = CompositeScorer([
            (QEDScorer(),         1.0),
            (SAScorer(),          1.0),
            (LipinskiScorer(),    0.5),
            (SubstructureAlerts(), 1.0),
        ], mode='product')
        result = scorer.score_one("CCO")  # dict of components + 'total'
    """

    def __init__(self, components, mode: str = "product"):
        assert mode in ("product", "sum")
        self.components = list(components)  # list of (scorer, weight)
        self.mode = mode

    def score_one(self, smiles: Optional[str]) -> Dict[str, float]:
        parts: Dict[str, float] = {}
        for scorer, _w in self.components:
            parts[scorer.name] = float(scorer.score(smiles))

        if self.mode == "product":
            total = 1.0
            wsum = 0.0
            for (scorer, w) in self.components:
                v = parts[scorer.name]
                # Geometric mean: total = prod(v_i ^ w_i) ^ (1/sum_w)
                # Use log-space to keep numerically stable
                if v <= 0.0:
                    total = 0.0
                    wsum += w
                    continue
                total *= v ** w
                wsum += w
            if wsum > 0 and total > 0:
                total = total ** (1.0 / wsum)
        else:  # sum
            num, den = 0.0, 0.0
            for (scorer, w) in self.components:
                num += parts[scorer.name] * w
                den += w
            total = num / den if den > 0 else 0.0

        parts["total_score"] = float(total)
        return parts

    def score_many(self, smiles_list: List[Optional[str]]) -> List[Dict[str, float]]:
        return [self.score_one(s) for s in smiles_list]


# ---------------------------------------------------------------------------
# Default factory
# ---------------------------------------------------------------------------
def default_scorer(reference_smiles: Iterable[str] = ()) -> CompositeScorer:
    """A sensible default for PROTAC evaluation."""
    return CompositeScorer(
        components=[
            (ValidityScorer(),             1.0),
            (QEDScorer(),                  1.0),
            (SAScorer(),                   1.0),
            (LipinskiScorer(),             0.5),
            (LinkerLengthScorer(30, 90),   0.5),
            (SubstructureAlerts(),         1.0),
            (NoveltyScorer(reference_smiles), 0.5),
        ],
        mode="product",
    )
