"""Project the 17 new ProtT5 embeddings into the existing protein UMAP.

The original UMAP reducer wasn't pickled, so a true transform() is not
available. Strategy: nearest-neighbor in 1024-d ProtT5 cosine space →
inherit the neighbor's 2D UMAP coordinates (with small jitter so points
don't perfectly stack) and cluster assignment.

Justification: UMAP preserves local cosine structure for its k-NN graph
(k=15 in the original fit). A protein that is the top cosine neighbor of
an existing protein lies inside the same local manifold patch, so the
neighbor's 2D coords are a faithful proxy for where UMAP would place it
in a refit. For only 17 added points among 37,391, a full refit would
shift the whole layout while preserving the same local outcome.

Output: protein_sequence_embedding.csv overwritten with
  - an added `is_new` column (0 for the original 37,391; 1 for new)
  - 17 new rows with the projected coords, NN-cluster, and paper metadata
The OLD file is preserved as protein_sequence_embedding.PRE_NEW.csv.
"""
from __future__ import annotations

import json
import re
from pathlib import Path

import h5py
import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent
CSV_IN = ROOT / "protein_sequence_embedding.csv"
CSV_OUT = CSV_IN
CSV_BACKUP = ROOT / "protein_sequence_embedding.PRE_NEW.csv"
H5_OLD = Path("/home/adsiordia/steroid_core_classifier/embeddings/all_steroid_uniprot_comprehensive_v2.h5")
H5_NEW = HERE / "new_proteins_18.h5"
NEW_FASTA = HERE / "new_proteins_18.clean.fasta"
NEW_FASTA_RAW = HERE / "new_proteins_18.fasta"

# rough mapping accession -> source paper, hand-built from new_proteins_v2.json
# and new_proteins_2024_2026.corrected.json. Used for tooltip / Annotation.
SOURCE_META = {
    "Q5LF84":     {"gene": "bsh (cbh)",        "paper": "Rimal et al., Nature 2024",
                   "doi": "10.1038/s41586-023-06990-w",
                   "function": "Bile salt hydrolase + amine acyltransferase"},
    "C8WL28":     {"gene": "Elen_2451 (fdhD)", "paper": "McCurry et al., Cell 2024",
                   "doi": "10.1016/j.cell.2024.05.005",
                   "function": "FdhD accessory of 21-dehydroxylase complex"},
    "C8WL29":     {"gene": "Elen_2452",        "paper": "McCurry et al., Cell 2024",
                   "doi": "10.1016/j.cell.2024.05.005",
                   "function": "4Fe-4S ferredoxin in 21-dehydroxylase"},
    "C8WL30":     {"gene": "Elen_2453",        "paper": "McCurry et al., Cell 2024",
                   "doi": "10.1016/j.cell.2024.05.005",
                   "function": "Mo-dependent 21-dehydroxylase catalytic subunit"},
    "C8WL31":     {"gene": "Elen_2454",        "paper": "McCurry et al., Cell 2024",
                   "doi": "10.1016/j.cell.2024.05.005",
                   "function": "SPFH/band-7 scaffold of 21-dehydroxylase"},
    "MFU7516964_1": {"gene": "OsrA", "paper": "Jacoby et al., Cell Host Microbe 2024",
                     "doi": "10.1016/j.chom.2025.09.014",
                     "function": "3-oxo-Delta1 reductase (prednisolone reduction)"},
    "MFU7517346_1": {"gene": "OsrB", "paper": "Jacoby et al., Cell Host Microbe 2024",
                     "doi": "10.1016/j.chom.2025.09.014",
                     "function": "3-oxo-Delta4 reductase (cortisol -> 5beta-DH)"},
    "MFU7515415_1": {"gene": "OsrC", "paper": "Jacoby et al., Cell Host Microbe 2024",
                     "doi": "10.1016/j.chom.2025.09.014",
                     "function": "3-oxo-5beta oxidoreductase"},
    "A0A8F5DVT9": {"gene": "bsh",      "paper": "Christensenella minuta BSH (multi-paper)",
                   "doi": "(multiple)",
                   "function": "Generic bile-salt hydrolase reference"},
    "B4YSU0":     {"gene": "baiCD",    "paper": "Cai et al., Nat Rev Gastro Hep 2024",
                   "doi": "10.1038/s41575-024-00914-3",
                   "function": "Bile-acid 7alpha-dehydroxylation enzyme"},
    "Q5LG86":     {"gene": "BF9343_1074", "paper": "Pathogenic B. fragilis (bioRxiv 2024)",
                   "doi": "10.1101/2024.06.19.599758",
                   "function": "B. fragilis BSH-cluster neighbor (uncharacterized)"},
    "Q5LF82":     {"gene": "BF9343_1435", "paper": "Pathogenic B. fragilis (bioRxiv 2024)",
                   "doi": "10.1101/2024.06.19.599758",
                   "function": "B. fragilis BSH-cluster neighbor (hypothetical exported)"},
    "Q5LF81":     {"gene": "BF9343_1436", "paper": "Pathogenic B. fragilis (bioRxiv 2024)",
                   "doi": "10.1101/2024.06.19.599758",
                   "function": "B. fragilis outer-membrane protein"},
    "A0A380Z384": {"gene": "BF9343_1437", "paper": "Pathogenic B. fragilis (bioRxiv 2024)",
                   "doi": "10.1101/2024.06.19.599758",
                   "function": "B. fragilis outer-membrane protein"},
    "A0A380YZ17": {"gene": "BF9343_1444 (mobC)", "paper": "Pathogenic B. fragilis (bioRxiv 2024)",
                   "doi": "10.1101/2024.06.19.599758",
                   "function": "Mobilisation protein C"},
    "Q5LDV5":     {"gene": "BF9343_1919", "paper": "Pathogenic B. fragilis (bioRxiv 2024)",
                   "doi": "10.1101/2024.06.19.599758",
                   "function": "B. fragilis uncharacterized"},
    "Q5LCD5":     {"gene": "BF9343_2449", "paper": "Pathogenic B. fragilis (bioRxiv 2024)",
                   "doi": "10.1101/2024.06.19.599758",
                   "function": "B. fragilis membrane protein"},
}


def read_fasta_records(path: Path) -> dict[str, tuple[str, str]]:
    """Returns acc -> (header_description, sequence)."""
    records: dict[str, tuple[str, str]] = {}
    acc = None
    desc = ""
    seq_parts: list[str] = []
    for line in path.read_text().splitlines():
        if line.startswith(">"):
            if acc:
                records[acc] = (desc, "".join(seq_parts))
            header = line.lstrip(">").strip()
            # strip the leading accession token, keep the rest as description
            toks = header.split(maxsplit=1)
            acc = toks[0].replace(".", "_").replace("/", "_")
            desc = toks[1] if len(toks) > 1 else ""
            seq_parts = []
        else:
            seq_parts.append(line.strip())
    if acc:
        records[acc] = (desc, "".join(seq_parts))
    return records


def main() -> None:
    print("Backing up existing CSV ...")
    df = pd.read_csv(CSV_IN)
    print(f"  current rows: {len(df):,}")
    df.to_csv(CSV_BACKUP, index=False)

    # --- 1. Load existing 1024-d embeddings, build acc-> vector map -----
    print("\nLoading existing h5 ...")
    with h5py.File(H5_OLD, "r") as f:
        old_keys = list(f.keys())
        old_mat = np.stack([f[k][:] for k in old_keys]).astype(np.float32)
    print(f"  shape: {old_mat.shape}")
    # h5 key format: "A0A016SKM3 A0A016SKM3_9BILA | Ancylostoma ceylanicum"
    # accession is the first whitespace-separated token
    old_acc = [k.split()[0] for k in old_keys]
    acc_to_idx = {a: i for i, a in enumerate(old_acc)}

    # --- 2. Load new embeddings -----------------------------------------
    print("\nLoading new h5 ...")
    with h5py.File(H5_NEW, "r") as f:
        new_keys = list(f.keys())
        new_mat = np.stack([f[k][:] for k in new_keys]).astype(np.float32)
    print(f"  shape: {new_mat.shape}")
    print(f"  accessions: {new_keys}")

    # --- 3. Normalize and find nearest neighbor in old set ---------------
    print("\nFinding nearest neighbors (cosine) ...")
    old_norm = old_mat / (np.linalg.norm(old_mat, axis=1, keepdims=True) + 1e-9)
    new_norm = new_mat / (np.linalg.norm(new_mat, axis=1, keepdims=True) + 1e-9)
    sims = new_norm @ old_norm.T                  # (17, 37391)
    nn_idx = np.argmax(sims, axis=1)
    nn_sim = sims[np.arange(len(new_keys)), nn_idx]

    print("\n  acc        -> nearest existing -> sim   -> cluster")
    for i, acc in enumerate(new_keys):
        nn_acc = old_acc[nn_idx[i]]
        nn_row = df.loc[df["Entry"] == nn_acc]
        cl = int(nn_row["clusters"].iloc[0]) if not nn_row.empty else -1
        print(f"  {acc:<13} -> {nn_acc:<12} -> {nn_sim[i]:.3f} -> cluster {cl}")

    # --- 4. Build new rows ---------------------------------------------
    print("\nBuilding new rows ...")
    fasta_records = read_fasta_records(NEW_FASTA)
    raw_records = read_fasta_records(NEW_FASTA_RAW)  # has full descriptions

    # build raw_acc -> raw header description map
    full_desc_map: dict[str, str] = {}
    with NEW_FASTA_RAW.open() as fh:
        for line in fh:
            if line.startswith(">"):
                s = line.lstrip(">").strip()
                # Normalize accession the same way the embedder did
                m = re.match(r"(?:sp|tr)\|([A-Z0-9]+)\|", s)
                if m:
                    acc = m.group(1)
                else:
                    acc = s.split()[0].replace(".", "_").replace("/", "_")
                full_desc_map[acc] = s

    rng = np.random.default_rng(seed=42)
    new_rows = []
    for i, acc in enumerate(new_keys):
        nn_row = df.loc[df["Entry"] == old_acc[nn_idx[i]]].iloc[0]
        meta = SOURCE_META.get(acc, {"gene": acc, "paper": "?", "doi": "?", "function": "?"})
        seq = fasta_records[acc][1]
        desc = full_desc_map.get(acc, fasta_records[acc][0])

        # organism: pull from FASTA header description (after OS=)
        m_os = re.search(r"OS=([^=]+?)(?:\s+OX=|\s+GN=|\s+PE=|\s+SV=|$)", desc)
        if m_os:
            organism = m_os.group(1).strip()
        else:
            # NCBI-style "[Clostridium sp. HCS.1]"
            m_brack = re.search(r"\[([^\]]+)\]", desc)
            organism = m_brack.group(1) if m_brack else ""

        # small Gaussian jitter so multiple new points sharing a NN don't stack
        jitter_x, jitter_y = rng.normal(0, 0.15, size=2)
        umap_1 = float(nn_row["UMAP_1"]) + jitter_x
        umap_2 = float(nn_row["UMAP_2"]) + jitter_y

        new_rows.append({
            "Unnamed: 0": len(df) + i,
            "Entry": acc.replace("_1", "") if acc.startswith("MFU") else acc,
            "Entry Name": acc,
            "Protein names": meta["function"],
            "Gene Names": meta["gene"],
            "Organism": organism,
            "Length": len(seq),
            "Sequence": seq,
            "Annotation": f'{meta["paper"]} | DOI: {meta["doi"]}',
            "source": "NEW_2024_2026",
            "ChEBI ID": "",
            "Rhea ID": "",
            "SMILES": "",
            "UMAP_1": umap_1,
            "UMAP_2": umap_2,
            "clusters": int(nn_row["clusters"]),
            "Compound Name": "",
        })

    # --- 5. Mark is_new on the combined frame ---------------------------
    df["is_new"] = 0
    new_df = pd.DataFrame(new_rows)
    new_df["is_new"] = 1
    combined = pd.concat([df, new_df], ignore_index=True)
    combined.to_csv(CSV_OUT, index=False)
    print(f"\nWrote {CSV_OUT}")
    print(f"  total rows: {len(combined):,}  ({int(combined['is_new'].sum())} new)")
    print(f"  backup at:  {CSV_BACKUP}")


if __name__ == "__main__":
    main()
