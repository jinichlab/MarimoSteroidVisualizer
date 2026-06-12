"""Refine the catalog triage.

Step 1 — for every candidate accession from a HAS_ACCESSIONS paper, hit
the live UniProt / NCBI service and only keep those that actually resolve
to a protein. Drops the regex's false positives (mouse strain B6C3F1,
buffer K2HPO4, PMC IDs caught by the locus pattern, etc.).

Step 2 — for each surviving candidate, fetch the FASTA and write it to
literature/sequences_v2/<safe_id>.fasta, plus add a structured record to
new_proteins_v2.json with explicit provenance:

    {
      "gene_name_or_label": ...,
      "organism_from_fasta": ...,
      "accession": ...,
      "accession_type": "uniprot" | "ncbi_protein" | "locus_tag",
      "validated": true,
      "fasta_path": "...",
      "length": ...,
      "paper_doi": ...,
      "paper_title": ...,
      "paper_text_source": "PMC:...",
      "evidence": "literal accession appeared in PMC full text"
    }

Step 3 — build manual_review_queue.tsv from the 28 ENZYME_LIKELY rows
plus the 5 catalog-flagged 'unknown' rows, with enzyme keyword summary so
a human can prioritize.
"""
from __future__ import annotations

import json
import re
import time
import urllib.parse
import urllib.request
from pathlib import Path

HERE = Path(__file__).resolve().parent
TRIAGE = HERE / "triage_catalog_report.json"
CATALOG = HERE / "dorrestein_catalog_2024_2026.json"
OUT_V2_JSON = HERE / "new_proteins_v2.json"
OUT_FASTA_DIR = HERE / "sequences_v2"
OUT_REVIEW_TSV = HERE / "manual_review_queue.tsv"
OUT_REJECT_TSV = HERE / "rejected_candidates.tsv"

OUT_FASTA_DIR.mkdir(exist_ok=True)

UA = "literature-refine (adsiordia@ucsd.edu)"
DELAY = 0.4

# Strings that pass the regex but are NOT proteins. Add as discovered.
STOP_STRINGS = {
    "K2HPO4", "Na2HPO4", "NaH2PO4", "MgSO4", "CaCl2", "KH2PO4",  # buffers
    "B6C3F1", "C57BL", "DBA2J", "BALB", "FVB1J", "129SVE",       # mouse strains
    "ATP", "NADP", "NADH", "ADP", "FADH", "FMNH",                # cofactors
    "PMC", "DOI", "PMID",
}
# PMC ID pattern caught by NCBI regex — kill them all
PMC_RX = re.compile(r"^PMC\d+\.\d+$", re.I)


def http_get(url: str, timeout: int = 30) -> tuple[str, int]:
    req = urllib.request.Request(url, headers={"User-Agent": UA})
    try:
        with urllib.request.urlopen(req, timeout=timeout) as r:
            return r.read().decode("utf-8", errors="replace"), r.status
    except urllib.error.HTTPError as e:
        return "", e.code
    except Exception:
        return "", -1


def parse_fasta(text: str) -> tuple[str, str]:
    text = text.strip()
    if not text.startswith(">"):
        return "", ""
    lines = text.split("\n")
    header = lines[0].lstrip(">").strip()
    seq = "".join(l.strip() for l in lines[1:] if not l.startswith(">"))
    seq = re.sub(r"[^A-Za-z*]", "", seq)
    return header, seq


def validate_uniprot(acc: str) -> tuple[bool, str, int, str]:
    """Returns (ok, header, length, raw_fasta)."""
    if acc in STOP_STRINGS:
        return False, "stop_string", 0, ""
    url = f"https://rest.uniprot.org/uniprotkb/{acc}.fasta"
    text, status = http_get(url)
    time.sleep(DELAY)
    if status != 200 or not text:
        return False, f"http_{status}", 0, ""
    header, seq = parse_fasta(text)
    # UniProt sometimes returns merge/deletion pages — make sure we got a real FASTA
    if not seq or len(seq) < 30:
        return False, "short_or_invalid", len(seq), text
    return True, header, len(seq), text


def validate_ncbi(acc: str) -> tuple[bool, str, int, str]:
    if acc in STOP_STRINGS or PMC_RX.match(acc):
        return False, "stop_string_or_pmc", 0, ""
    url = (
        "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi"
        f"?db=protein&id={urllib.parse.quote(acc)}&rettype=fasta&retmode=text"
    )
    text, status = http_get(url, timeout=45)
    time.sleep(DELAY)
    if status != 200 or not text:
        return False, f"http_{status}", 0, ""
    header, seq = parse_fasta(text)
    if not seq or len(seq) < 30:
        return False, "short_or_invalid", len(seq), text
    return True, header, len(seq), text


def validate_locus_via_uniprot_search(locus: str, organism_hint: str = "") -> tuple[bool, str, int, str, str]:
    """Try to look up a locus tag via UniProt search; return (ok, header, length, raw_fasta, accession)."""
    if locus in STOP_STRINGS:
        return False, "stop", 0, "", ""
    q = f'gene:"{locus}"'
    if organism_hint:
        q += f' AND organism_name:"{organism_hint.split()[0]}"'
    url = (
        "https://rest.uniprot.org/uniprotkb/search?"
        f"query={urllib.parse.quote(q)}&format=fasta&size=1"
    )
    text, status = http_get(url, timeout=30)
    time.sleep(DELAY)
    if status != 200 or not text or not text.startswith(">"):
        return False, f"http_{status}", 0, "", ""
    header, seq = parse_fasta(text)
    if not seq or len(seq) < 30:
        return False, "short", len(seq), text, ""
    # extract accession from header like "sp|Q5LF84|..."
    m = re.search(r"\|([A-Z0-9]{6,10})\|", header)
    acc = m.group(1) if m else "?"
    return True, header, len(seq), text, acc


def safe_filename(s: str) -> str:
    return re.sub(r"[^A-Za-z0-9_.-]+", "_", s).strip("_")


def main() -> None:
    triage = json.loads(TRIAGE.read_text())
    catalog = {p["doi"]: p for p in json.loads(CATALOG.read_text())["papers"] if p.get("doi")}

    has_acc_papers = [r for r in triage if r["verdict"] == "HAS_ACCESSIONS"]
    enzyme_likely = [r for r in triage if r["verdict"] == "ENZYME_LIKELY"]

    print(f"HAS_ACCESSIONS papers to validate: {len(has_acc_papers)}")
    print(f"ENZYME_LIKELY papers to queue for review: {len(enzyme_likely)}\n")

    validated_records = []
    rejected_records = []

    for i, paper in enumerate(has_acc_papers, 1):
        doi = paper["doi"]
        title = paper["title"]
        print(f"\n[{i}/{len(has_acc_papers)}] {title[:90]}")
        print(f"  DOI: {doi}")

        # Collect candidates
        cand_set = []
        for raw, kind in [
            (paper.get("uniprot_candidates", ""), "uniprot"),
            (paper.get("ncbi_candidates", ""), "ncbi_protein"),
            (paper.get("locus_candidates", ""), "locus_tag"),
        ]:
            for c in raw.split(";"):
                c = c.strip()
                if c and c not in STOP_STRINGS:
                    cand_set.append((c, kind))

        # Dedup, preserve order
        seen = set()
        cands = []
        for c, k in cand_set:
            if c not in seen:
                seen.add(c)
                cands.append((c, k))

        if not cands:
            continue

        # organism hint: pull from title — often something like "...in Eggerthella..." or "...Bacteroides fragilis..."
        organism_hint = ""
        m_org = re.search(
            r"\b(Bacteroides|Eggerthella|Bifidobacterium|Clostridium|Lactobacillus|Enterococcus|Escherichia|Mycobacterium|Streptococcus|Akkermansia|Faecalibacterium|Lactococcus|Ruminococcus|Anaerostipes|Blautia|Roseburia)\s+\w+",
            title,
            re.I,
        )
        if m_org:
            organism_hint = m_org.group(0)

        for c, kind in cands:
            print(f"   validating {kind}:{c} ...", end=" ", flush=True)
            if kind == "uniprot":
                ok, header, length, fasta = validate_uniprot(c)
                final_acc = c
            elif kind == "ncbi_protein":
                ok, header, length, fasta = validate_ncbi(c)
                final_acc = c
            else:  # locus_tag
                ok, header, length, fasta, mapped_acc = validate_locus_via_uniprot_search(c, organism_hint)
                final_acc = mapped_acc if ok else c

            if ok:
                print(f"OK len={length}")
                fname = safe_filename(f"{c}_{final_acc}") + ".fasta"
                (OUT_FASTA_DIR / fname).write_text(fasta if fasta.endswith("\n") else fasta + "\n")
                validated_records.append({
                    "candidate": c,
                    "accession_type": kind,
                    "validated_accession": final_acc,
                    "fasta_header": header[:200],
                    "length": length,
                    "fasta_path": f"sequences_v2/{fname}",
                    "paper_doi": doi,
                    "paper_title": title,
                    "paper_text_source": paper.get("text_source", ""),
                    "evidence": (
                        f"Candidate '{c}' surfaced from {paper.get('text_source','?')} for paper {doi}; "
                        f"validated by live fetch from "
                        f"{'UniProt' if kind in ('uniprot','locus_tag') else 'NCBI'} "
                        f"resolving to {final_acc} ({length} aa)."
                    ),
                })
            else:
                print(f"REJECT ({header})")
                rejected_records.append({
                    "candidate": c,
                    "kind": kind,
                    "reason": header,
                    "paper_doi": doi,
                    "paper_title": title,
                })

    # Manual review queue: ENZYME_LIKELY + catalog 'unknown'
    catalog_unknowns = [
        {
            "doi": p.get("doi", ""),
            "year": p.get("year", ""),
            "journal": p.get("journal_or_preprint_server", ""),
            "title": p.get("title", ""),
            "verdict": "catalog_unknown",
            "enzyme_keyword_hits": "",
        }
        for p in json.loads(CATALOG.read_text())["papers"]
        if str(p.get("introduces_novel_protein", "")).lower() == "unknown"
    ]
    queue = []
    for r in enzyme_likely + catalog_unknowns:
        queue.append({
            "doi": r.get("doi", ""),
            "year": r.get("year", ""),
            "journal": r.get("journal", ""),
            "title": r.get("title", ""),
            "verdict": r.get("verdict", ""),
            "enzyme_keyword_hits": r.get("enzyme_keyword_hits", ""),
        })

    OUT_V2_JSON.write_text(json.dumps(validated_records, indent=2))
    cols = ["doi", "year", "journal", "verdict", "enzyme_keyword_hits", "title"]
    with OUT_REVIEW_TSV.open("w") as fh:
        fh.write("\t".join(cols) + "\n")
        for r in queue:
            fh.write("\t".join(str(r.get(c, "")).replace("\t", " ") for c in cols) + "\n")
    with OUT_REJECT_TSV.open("w") as fh:
        fh.write("candidate\tkind\treason\tpaper_doi\tpaper_title\n")
        for r in rejected_records:
            fh.write("\t".join(str(r.get(c, "")).replace("\t", " ") for c in
                               ["candidate", "kind", "reason", "paper_doi", "paper_title"]) + "\n")

    print("\n\n=== DONE ===")
    print(f"Validated proteins: {len(validated_records)}")
    print(f"Rejected candidates: {len(rejected_records)}")
    print(f"Manual review queue: {len(queue)}")
    print(f"\nOutputs:")
    print(f"  {OUT_V2_JSON}")
    print(f"  {OUT_FASTA_DIR}/  ({len(validated_records)} FASTAs)")
    print(f"  {OUT_REVIEW_TSV}")
    print(f"  {OUT_REJECT_TSV}")


if __name__ == "__main__":
    main()
