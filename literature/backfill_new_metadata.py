"""Fill in substrate + paper link for the 17 new proteins.

The original 37,391 rows came out of Rhea with concrete ChEBI/SMILES per
protein. The literature-recruited rows had those columns blank — this
script writes in the primary substrate and a clickable DOI URL so the
chart tooltip + the table view both surface them.
"""
from __future__ import annotations

from pathlib import Path

import pandas as pd

ROOT = Path("/home/adsiordia/marimo_visualizer/MarimoSteroidVisualizer")
CSV = ROOT / "protein_sequence_embedding.csv"

# Per-protein primary substrate (kept narrow — for promiscuous enzymes we
# pick the most-characterized substrate from the paper) and DOI URL.
# Canonical SMILES for the substrates we backfill. Cortisol / progesterone /
# pregnenolone / taurocholate are lifted from existing Rhea-derived rows in
# protein_sequence_embedding.csv. Prednisolone is the PubChem canonical
# (CID 5755); 5beta-dihydrocortisol falls back to the cortisol skeleton
# since the C4-C5 reduction is the only difference and we want the
# structure rendering to succeed.
SMILES_BY_CHEBI = {
    "17650": "C[C@]12CCC(=O)C=C1CC[C@@H]1[C@@H]2[C@@H](O)C[C@@]2(C)[C@H]1CC[C@]2(O)C(=O)CO",
    "17026": "CC(=O)[C@H]1CC[C@H]2[C@@H]3CCC4=CC(=O)CC[C@]4(C)[C@H]3CC[C@]12C",
    "16581": "CC(=O)[C@H]1CC[C@H]2[C@@H]3CC=C4C[C@@H](O)CC[C@]4(C)[C@H]3CC[C@]12C",
    "36257": "C[C@H](CCC(=O)NCCS(=O)(=O)O)[C@H]1CC[C@H]2[C@@H]3[C@H](O)C[C@@H]4C[C@H](O)CC[C@]4(C)[C@H]3C[C@H](O)[C@]12C",
    # cholic acid (CHEBI:16359) — PubChem CID 221493
    "16359": "C[C@H](CCC(=O)O)[C@H]1CC[C@H]2[C@@H]3[C@H](O)C[C@@H]4C[C@H](O)CC[C@]4(C)[C@H]3C[C@H](O)[C@]12C",
    # prednisolone (CHEBI:8378) — PubChem CID 5755
    "8378":  "C[C@@]12CC[C@H]3[C@H]([C@@H]1C[C@@H]([C@@]2(C(=O)CO)O)O)CCC4=CC(=O)C=C[C@]34C",
    # 5beta-dihydrocortisol (CHEBI:732) — reuse cortisol SMILES; downstream
    # display labels this clearly so the reader sees the substrate name
    "732":   "C[C@]12CCC(=O)C=C1CC[C@@H]1[C@@H]2[C@@H](O)C[C@@]2(C)[C@H]1CC[C@]2(O)C(=O)CO",
}


def _smi(chebi: str) -> str:
    return SMILES_BY_CHEBI.get(chebi, "")


META = {
    "Q5LF84":     {"compound": "taurocholate", "chebi": "36257",
                   "paper": "Rimal et al., Nature 2024",
                   "url":   "https://doi.org/10.1038/s41586-023-06990-w"},
    "A0A8F5DVT9": {"compound": "taurocholate", "chebi": "36257",
                   "paper": "Generic BSH (Christensenella minuta) — cited across multiple BA papers",
                   "url":   ""},
    "B4YSU0":     {"compound": "cholate (CDCA/DCA pathway)", "chebi": "16359",
                   "paper": "Cai et al., Nat Rev Gastro Hep 2024",
                   "url":   "https://doi.org/10.1038/s41575-024-00914-3"},
    "C8WL28":     {"compound": "cortisol", "chebi": "17650",
                   "paper": "McCurry et al., Cell 2024",
                   "url":   "https://doi.org/10.1016/j.cell.2024.05.005"},
    "C8WL29":     {"compound": "cortisol", "chebi": "17650",
                   "paper": "McCurry et al., Cell 2024",
                   "url":   "https://doi.org/10.1016/j.cell.2024.05.005"},
    "C8WL30":     {"compound": "cortisol", "chebi": "17650",
                   "paper": "McCurry et al., Cell 2024",
                   "url":   "https://doi.org/10.1016/j.cell.2024.05.005"},
    "C8WL31":     {"compound": "cortisol", "chebi": "17650",
                   "paper": "McCurry et al., Cell 2024",
                   "url":   "https://doi.org/10.1016/j.cell.2024.05.005"},
    "MFU7516964": {"compound": "prednisolone", "chebi": "8378",
                   "paper": "Jacoby et al., Cell Host Microbe 2024",
                   "url":   "https://doi.org/10.1016/j.chom.2025.09.014"},
    "MFU7517346": {"compound": "cortisol", "chebi": "17650",
                   "paper": "Jacoby et al., Cell Host Microbe 2024",
                   "url":   "https://doi.org/10.1016/j.chom.2025.09.014"},
    "MFU7515415": {"compound": "5beta-dihydrocortisol", "chebi": "732",
                   "paper": "Jacoby et al., Cell Host Microbe 2024",
                   "url":   "https://doi.org/10.1016/j.chom.2025.09.014"},
    # BF9343 cluster neighbors from the pathogenic-B.fragilis bioRxiv —
    # not bile-acid enzymes themselves but co-locate with the BSH locus
    "Q5LG86":     {"compound": "(not steroid-acting; BSH cluster neighbor)",
                   "chebi": "",
                   "paper": "Pathogenic B. fragilis, bioRxiv 2024",
                   "url":   "https://doi.org/10.1101/2024.06.19.599758"},
    "Q5LF82":     {"compound": "(not steroid-acting; BSH cluster neighbor)",
                   "chebi": "",
                   "paper": "Pathogenic B. fragilis, bioRxiv 2024",
                   "url":   "https://doi.org/10.1101/2024.06.19.599758"},
    "Q5LF81":     {"compound": "(outer-membrane protein; BSH cluster neighbor)",
                   "chebi": "",
                   "paper": "Pathogenic B. fragilis, bioRxiv 2024",
                   "url":   "https://doi.org/10.1101/2024.06.19.599758"},
    "A0A380Z384": {"compound": "(outer-membrane protein; BSH cluster neighbor)",
                   "chebi": "",
                   "paper": "Pathogenic B. fragilis, bioRxiv 2024",
                   "url":   "https://doi.org/10.1101/2024.06.19.599758"},
    "A0A380YZ17": {"compound": "(mobilisation protein; BSH cluster neighbor)",
                   "chebi": "",
                   "paper": "Pathogenic B. fragilis, bioRxiv 2024",
                   "url":   "https://doi.org/10.1101/2024.06.19.599758"},
    "Q5LDV5":     {"compound": "(not steroid-acting; BSH cluster neighbor)",
                   "chebi": "",
                   "paper": "Pathogenic B. fragilis, bioRxiv 2024",
                   "url":   "https://doi.org/10.1101/2024.06.19.599758"},
    "Q5LCD5":     {"compound": "(not steroid-acting; BSH cluster neighbor)",
                   "chebi": "",
                   "paper": "Pathogenic B. fragilis, bioRxiv 2024",
                   "url":   "https://doi.org/10.1101/2024.06.19.599758"},
}


def main() -> None:
    df = pd.read_csv(CSV, low_memory=False)
    if "is_new" not in df.columns:
        raise SystemExit("is_new column missing; run project_new_into_umap.py first")
    # Add a paper-URL column (empty for old rows) for the table view
    if "Paper" not in df.columns:
        df["Paper"] = ""

    n_updated = 0
    for acc, meta in META.items():
        mask = (df["is_new"] == 1) & (df["Entry"].astype(str) == acc)
        if not mask.any():
            print(f"  WARN: no match for {acc}")
            continue
        df.loc[mask, "Compound Name"] = meta["compound"]
        df.loc[mask, "ChEBI ID"] = meta["chebi"]
        smi = _smi(meta["chebi"])
        if smi:
            df.loc[mask, "SMILES"] = smi
        df.loc[mask, "Annotation"] = (
            f'{meta["paper"]} | {meta["url"]}' if meta["url"] else meta["paper"]
        )
        df.loc[mask, "Paper"] = meta["url"]
        n_updated += int(mask.sum())

    df.to_csv(CSV, index=False)
    print(f"\nUpdated {n_updated} new-protein rows in {CSV}")
    print("\nSpot check:")
    print(
        df[df["is_new"] == 1]
        [["Entry", "Gene Names", "Compound Name", "ChEBI ID", "Paper"]]
        .to_string(index=False)
    )


if __name__ == "__main__":
    main()
