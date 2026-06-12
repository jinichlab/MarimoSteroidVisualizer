"""Render new_molecules_2024_2026.json as labeled 2D grids for human review.

Outputs one PNG per molecule class plus a combined PNG. Run from anywhere:
    LD_LIBRARY_PATH=/home/adsiordia/miniconda3/lib \
    /home/adsiordia/miniconda3/bin/python \
    /home/adsiordia/marimo_visualizer/MarimoSteroidVisualizer/literature/render_molecules_grid.py
"""
import json
from collections import defaultdict
from pathlib import Path

from rdkit import Chem
from rdkit.Chem import AllChem, Draw

HERE = Path(__file__).resolve().parent
SRC = HERE / "new_molecules_2024_2026.json"
OUT_DIR = HERE / "molecule_grids"
OUT_DIR.mkdir(exist_ok=True)

with SRC.open() as fh:
    mols_json = json.load(fh)

by_class: dict[str, list[dict]] = defaultdict(list)
for entry in mols_json:
    by_class[entry["class"]].append(entry)

def parse(entry):
    mol = Chem.MolFromSmiles(entry["smiles"])
    if mol is None:
        return None
    AllChem.Compute2DCoords(mol)
    src_tag = "PubChem" if entry["smiles_source"] == "pubchem" else "constructed"
    name = entry["name"]
    if len(name) > 38:
        name = name[:35] + "..."
    mol.SetProp("_Name", f"{name}\n[{src_tag}]")
    return mol

all_mols, all_labels = [], []
for klass, entries in by_class.items():
    parsed = [(e, parse(e)) for e in entries]
    bad = [e["name"] for e, m in parsed if m is None]
    good = [(e, m) for e, m in parsed if m is not None]
    if bad:
        print(f"[{klass}] SMILES failed to parse: {bad}")
    if not good:
        continue
    mols = [m for _, m in good]
    labels = [m.GetProp("_Name") for m in mols]
    img = Draw.MolsToGridImage(
        mols,
        molsPerRow=3,
        subImgSize=(360, 320),
        legends=labels,
        useSVG=False,
    )
    out = OUT_DIR / f"{klass.replace(' ', '_').replace('/', '_')}.png"
    img.save(out)
    print(f"[{klass}] {len(mols)} molecules -> {out}")
    all_mols.extend(mols)
    all_labels.extend([f"{e['class']}: {lab}" for (e, _), lab in zip(good, labels)])

combined = Draw.MolsToGridImage(
    all_mols,
    molsPerRow=3,
    subImgSize=(360, 320),
    legends=all_labels,
    useSVG=False,
)
combined_out = OUT_DIR / "ALL_new_molecules.png"
combined.save(combined_out)
print(f"\nCombined ({len(all_mols)} molecules) -> {combined_out}")
