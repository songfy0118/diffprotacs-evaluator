"""
make_demo_data.py
=================
Generate a small fake DiffPROTACs output directory by taking known
SMILES, embedding them in 3D with RDKit, and saving as .xyz in the
same layout DiffPROTACs writes (<output_dir>/<uuid>/<i>_.xyz).

This lets you sanity-check the evaluator without running the diffusion
model. Run:
    python make_demo_data.py
to create ./demo_output with ~8 molecules.
"""

import os
import shutil
from pathlib import Path

from rdkit import Chem
from rdkit.Chem import AllChem


# A small mix: a couple of real PROTAC-style molecules, some normal drugs,
# and an intentionally weird/invalid one.
DEMO_SMILES = [
    # Neutral, no carboxylic acids / phenols (which xyz2mol struggles
    # with at total_charge=0 because the H is implicit). All of these
    # parse cleanly back from xyz in a neutral Lewis structure.
    ("benzene",       "c1ccccc1"),
    ("toluene",       "Cc1ccccc1"),
    ("naphthalene",   "c1ccc2ccccc2c1"),
    ("anisole",       "COc1ccccc1"),
    ("caffeine",      "Cn1cnc2c1c(=O)n(C)c(=O)n2C"),
    # Small PROTAC-like fragment chains with amide linkers (neutral)
    ("amide_linker1", "CC(=O)NCCCCCNC(=O)C"),
    ("amide_linker2", "Cc1ccc(NC(=O)CCNC(=O)C)cc1"),
    ("peg_like",      "COCCOCCOCCOC"),
]


def smi_to_xyz_block(smi: str) -> str:
    """SMILES -> 3D coordinates -> xyz file content (DiffPROTACs format).

    Note: We do NOT aggressively optimize geometry, because too-tight
    optimization can collapse aromatic bond lengths and break xyz2mol
    bond inference downstream. ETKDG already gives chemically sensible
    starting geometries.
    """
    mol = Chem.MolFromSmiles(smi)
    if mol is None:
        return ""
    mol = Chem.AddHs(mol)
    params = AllChem.ETKDGv3()
    params.randomSeed = 42
    if AllChem.EmbedMolecule(mol, params) != 0:
        return ""
    # Light UFF cleanup only — preserves bond lengths better than MMFF
    # for arbitrary druglike molecules.
    try:
        AllChem.UFFOptimizeMolecule(mol, maxIters=50)
    except Exception:
        pass
    mol = Chem.RemoveHs(mol)  # DiffPROTACs writes heavy atoms only
    conf = mol.GetConformer()
    lines = [f"{mol.GetNumAtoms()}", ""]  # count + blank comment
    for atom in mol.GetAtoms():
        pos = conf.GetAtomPosition(atom.GetIdx())
        lines.append(f"{atom.GetSymbol()} {pos.x:.6f} {pos.y:.6f} {pos.z:.6f}")
    return "\n".join(lines) + "\n"


def main(out_dir: str = "demo_output"):
    out = Path(out_dir)
    if out.exists():
        shutil.rmtree(out)
    out.mkdir(parents=True)

    n_written = 0
    for uuid, smi in DEMO_SMILES:
        mol_dir = out / uuid
        mol_dir.mkdir()
        # Two samples per uuid, mimicking DiffPROTACs' n_stability_samples
        for i in range(2):
            block = smi_to_xyz_block(smi)
            if not block:
                continue
            with open(mol_dir / f"{i}_.xyz", "w") as f:
                f.write(block)
            n_written += 1
        # Also write a true_ reference (will be skipped by collector)
        block = smi_to_xyz_block(smi)
        if block:
            with open(mol_dir / "true_.xyz", "w") as f:
                f.write(block)

    print(f"wrote {n_written} sample .xyz files under {out}/")
    print("structure:")
    for p in sorted(out.rglob("*.xyz")):
        print(f"  {p}")


if __name__ == "__main__":
    # Always write under the script's own directory, regardless of cwd
    here = os.path.dirname(os.path.abspath(__file__))
    main(os.path.join(here, "demo_output"))
