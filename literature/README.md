# Literature pipeline (2024-2026 recruitment)

This folder is the **second iteration** of the steroid visualizer dataset.
The original visualizer (`smiles_visualizer.py` at the repo root) was
built from a Rhea snapshot containing 37,391 proteins. This pipeline
recruits **17 additional proteins** characterized in 2024-2026 literature
that were not in that Rhea snapshot, projects them into the existing
UMAP, and renders them as star markers in the visualizer.

## What lives here

| Stage | Script | Output |
|---|---|---|
| 1. Catalog discovery | (manual) | `dorrestein_catalog_2024_2026.json` — 144 candidate papers |
| 2. Triage all papers | `triage_all_catalog.py` | `triage_catalog_report.json/tsv` — verdicts by signal strength |
| 3. Validate accessions | `refine_triage.py` | `new_proteins_v2.json`, `sequences_v2/*.fasta`, `manual_review_queue.tsv`, `rejected_candidates.tsv` |
| 4. Provenance audit | `audit_protein_provenance.py` | `new_proteins_2024_2026.audited.json`, `provenance_report.tsv` |
| 5. Fetch round-1 FASTAs | `fetch_protein_sequences.py` | `sequences/*.fasta`, `sequence_fetch_report.tsv` |
| 6. Build joint FASTA | `normalize_fasta.py` | `new_proteins_18.fasta` → `new_proteins_18.clean.fasta` |
| 7. ProtT5 embed | (external script) | `new_proteins_18.h5` (17 × 1024) |
| 8. Project into UMAP | `project_new_into_umap.py` | overwrites `../protein_sequence_embedding.csv` with 17 new rows tagged `is_new=1` |
| 9. Backfill substrate + paper | `backfill_new_metadata.py` | fills `Compound Name`, `ChEBI ID`, `SMILES`, `Annotation`, `Paper` for the new rows |

## Final dataset

18 candidates were considered. After dedup (Q5LF84 appears in both the
Rimal *Nature* paper and the *Pathogenic B. fragilis* bioRxiv as the same
gene), **17 unique proteins** were embedded:

- **8 from manual extraction (round 1)** — `new_proteins_2024_2026.json`
  + `new_proteins_2024_2026.corrected.json` (note: the McCurry DOI was
  initially wrong; corrected file is authoritative). Sources: Rimal et al.
  Nature 2024, McCurry et al. Cell 2024, Jacoby et al. Cell Host Microbe
  2024/25.
- **9 net new from the full-catalog sweep (round 2)** — `new_proteins_v2.json`,
  after dropping 6 regex-true-but-biology-false matches (e.g. pig trypsin
  from a mass-spec methods paper) and 1 duplicate (Q5LF84). Sources:
  *Pathogenic B. fragilis* (bioRxiv), Cai et al. *Nat Rev Gastro Hep* 2024,
  multi-paper *Christensenella minuta* BSH.
- **P32370** (baiH, C. scindens) was the only catalog hit already present
  in the original 37,391-protein dataset — it is **not** re-added.

A breakdown of biological relevance: 10 of the 17 are confirmed
steroid-acting enzymes (BSHs, baiCD, baiH paralog, OsrABC,
Elen_2451-2454); the remaining 7 are *Bacteroides fragilis* NCTC 9343
BSH-cluster neighbors (OMPs, mobC, conserved hypotheticals) that
co-locate with the BSH gene but do not themselves catalyse steroid
chemistry. They are kept in the dataset for cluster-context reasons but
have empty `SMILES` / `ChEBI ID` so the substrate rendering correctly
shows "[Invalid SMILES]" for them.

## Reproducing

The pipeline needs:

1. **ProtT5 model** (`Rostlab/prot_t5_xl_half_uniref50-enc`) accessible
   via HuggingFace. Use the dedicated venv to avoid the transformers 5.x
   SentencePiece breakage:
   ```bash
   /home/adsiordia/miniconda3/bin/python -m venv /path/to/prott5_env
   /path/to/prott5_env/bin/pip install \
       'transformers==4.46.3' 'tokenizers<0.21' \
       sentencepiece protobuf h5py torch
   ```
2. **The original ProtT5 h5** of the 37,391-protein Rhea set, located on
   the GPU server at
   `/home/adsiordia/steroid_core_classifier/embeddings/all_steroid_uniprot_comprehensive_v2.h5`.
   This file is ~150 MB and is **not** tracked in git — see
   `../UMAP data/` for the path used by the projection script. Ask
   Andre / the lab for a copy if reproducing from scratch.

End-to-end run:
```bash
cd literature/
# Re-validate triage candidates (live UniProt/NCBI fetch)
python refine_triage.py

# Normalize FASTA headers, embed via ProtT5
python normalize_fasta.py
LD_LIBRARY_PATH=/home/adsiordia/miniconda3/lib \
  /path/to/prott5_env/bin/python /path/to/prott5_embedder.py \
    -i new_proteins_18.clean.fasta \
    -o new_proteins_18.h5 --per_protein 1

# Project into existing UMAP and backfill substrate metadata
python project_new_into_umap.py
python backfill_new_metadata.py
```

After step 3, restart the visualizer:
```bash
cd .. && bash run_visualizer.sh
```

## How the new proteins surface in the visualizer

The visualizer (`../smiles_visualizer.py`) reads `is_new` from the
protein CSV and:
- Renders new proteins as **5-pointed stars** instead of circles, at
  ~12× the area with a black stroke.
- Surfaces a **gold-bordered "★ Newly recruited" callout** in the
  per-protein detail panel, with the source citation inline.
- Adds a **clickable "Source paper" DOI link** in the same panel.
- Includes `Entry` (the UniProt accession) in the search box, so users
  can type e.g. `Q5LF84` to highlight a specific star.

## Known caveats

- The new proteins are placed via **1-NN in 1024-d ProtT5 cosine space**,
  inheriting the neighbor's 2D UMAP coords with a small jitter. This is
  a good local approximation but is *not* equivalent to refitting the
  full UMAP. If you need fully-comparable coordinates, refit the whole
  thing using `../build_protein_centric_v2.py` after appending the new
  embeddings to the source h5.
- `MFU7515415` (OsrC, 5β-dihydrocortisol substrate) reuses the cortisol
  SMILES for rendering since the 5β-reduction adds two H's that don't
  change the macro structure visibly at the displayed scale. The label
  still says "5beta-dihydrocortisol" so the reader isn't misled.
- The "generic BSH" entry `A0A8F5DVT9` (Christensenella minuta) was
  resolved 7× during catalog triage because 7 different papers mentioned
  the term "bsh" generically. The deduped record is intentional — there
  is only one ProtT5 embedding for it.
