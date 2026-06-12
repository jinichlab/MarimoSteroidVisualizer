"""Audit protein-ID provenance for new_proteins_2024_2026.json.

For each protein entry: fetch the paper's PMC/Nature URL, search the full
text for the claimed UniProt accession AND NCBI/locus identifier, and
classify provenance into one of four buckets:

  IN_PAPER       - the literal accession appears in the fetched paper text
  INFERRED       - accession NOT in paper; was looked up by gene + organism on UniProt
  GUESS          - JSON itself labels it "community ortholog" or family-level
  UNVERIFIABLE   - paper text could not be fetched (paywall / non-PMC)

Also greps for nearby accession-like strings around the gene name, so we
can spot cases where the paper DID cite a UniProt ID that the extractor
missed.

Writes:
  new_proteins_2024_2026.audited.json   - original entries + 'provenance' fields
  provenance_report.tsv                  - human-scannable summary
"""
from __future__ import annotations

import json
import re
import sys
import time
import urllib.request
import urllib.parse
from pathlib import Path

HERE = Path(__file__).resolve().parent
SRC = HERE / "new_proteins_2024_2026.json"
OUT_JSON = HERE / "new_proteins_2024_2026.audited.json"
OUT_TSV = HERE / "provenance_report.tsv"

UA = "literature-audit (adsiordia@ucsd.edu)"
FETCH_DELAY = 0.6  # be polite to PMC / publisher servers

# Patterns for finding accession-like strings in arbitrary text.
UNIPROT_RX = re.compile(
    r"\b("
    r"[OPQ][0-9][A-Z0-9]{3}[0-9]"                       # Q5LF84, P12345
    r"|[A-NR-Z][0-9](?:[A-Z][A-Z0-9]{2}[0-9]){1,2}"     # C8WL28, A0A024R161
    r")\b"
)
NCBI_PROT_RX = re.compile(
    r"\b("
    r"(?:WP_|NP_|YP_|XP_|MFU|MFV)\d+(?:\.\d+)?"           # WP_013410903.1, MFU7516964.1
    r"|[A-Z]{3}\d{5,}\.\d+"                                # ACV56407.1
    r")\b"
)


def fetch_url(url: str, timeout: int = 30) -> tuple[str, int]:
    """Return (text, status). Empty text means the fetch failed."""
    req = urllib.request.Request(url, headers={"User-Agent": UA})
    try:
        with urllib.request.urlopen(req, timeout=timeout) as r:
            return r.read().decode("utf-8", errors="replace"), r.status
    except Exception as e:
        return "", -1


def strip_html(html: str) -> str:
    """Crude tag-strip — keeps a lot of whitespace/noise but is enough for grep."""
    html = re.sub(r"<script.*?</script>", " ", html, flags=re.S | re.I)
    html = re.sub(r"<style.*?</style>", " ", html, flags=re.S | re.I)
    html = re.sub(r"<[^>]+>", " ", html)
    # html entities we care about
    html = html.replace("&nbsp;", " ").replace("&amp;", "&").replace("&gt;", ">").replace("&lt;", "<")
    return re.sub(r"\s+", " ", html)


def extract_bare_accession(field: str) -> str | None:
    """Pull the first UniProt or NCBI-shaped accession out of a JSON field
    that may carry parenthetical notes, e.g. 'Q5DWB7 (community ortholog; ...)'."""
    if not field:
        return None
    m = UNIPROT_RX.search(field)
    if m:
        return m.group(1)
    m = NCBI_PROT_RX.search(field)
    if m:
        return m.group(1)
    return None


def context_window(text: str, needle: str, window: int = 200) -> str:
    """Return a snippet of text around the first occurrence of needle."""
    idx = text.lower().find(needle.lower())
    if idx < 0:
        return ""
    start = max(0, idx - window)
    end = min(len(text), idx + len(needle) + window)
    return text[start:end].strip()


def near_gene_accessions(text: str, gene_token: str, window: int = 400) -> list[str]:
    """Find UniProt/NCBI-style accessions within `window` chars of every
    occurrence of the gene name. These are candidate IDs the original
    extractor might have missed."""
    found: set[str] = set()
    gene_l = gene_token.lower()
    text_l = text.lower()
    start = 0
    while True:
        idx = text_l.find(gene_l, start)
        if idx < 0:
            break
        snippet = text[max(0, idx - window): idx + len(gene_token) + window]
        for m in UNIPROT_RX.finditer(snippet):
            found.add(m.group(1))
        for m in NCBI_PROT_RX.finditer(snippet):
            found.add(m.group(1))
        start = idx + len(gene_token)
    return sorted(found)


def doi_to_pmc_url(doi: str) -> str | None:
    """Use NCBI idconv to translate a paper DOI -> PMC HTML URL.
    Returns None if no PMC version is indexed."""
    if not doi:
        return None
    api = (
        "https://www.ncbi.nlm.nih.gov/pmc/utils/idconv/v1.0/"
        f"?ids={urllib.parse.quote(doi)}&format=json&tool=audit&email=adsiordia@ucsd.edu"
    )
    text, status = fetch_url(api)
    if not text:
        return None
    try:
        payload = json.loads(text)
        for rec in payload.get("records", []):
            pmcid = rec.get("pmcid")
            if pmcid:
                return f"https://pmc.ncbi.nlm.nih.gov/articles/{pmcid}/"
    except Exception:
        pass
    return None


def resolve_paper_url(entry: dict) -> tuple[str, str]:
    """Pick the best URL to fetch the actual paper text. Prefers PMC.
    Returns (url, source_label)."""
    src = (entry.get("source_url") or "").strip()
    doi = (entry.get("paper_DOI") or "").strip()

    # If the source_url is a database page (uniprot.org, ncbi.nlm.nih.gov/protein),
    # it isn't the paper — we want to fetch the paper itself via PMC if possible.
    src_l = src.lower()
    is_paper_url = (
        "pmc.ncbi.nlm.nih.gov/articles/" in src_l
        or "doi.org/" in src_l
        or "nature.com/articles/" in src_l
        or "cell.com/" in src_l
        or "sciencedirect.com/" in src_l
    )
    if is_paper_url:
        return src, "source_url"

    # Try DOI -> PMC translation
    if doi:
        pmc_url = doi_to_pmc_url(doi)
        if pmc_url:
            return pmc_url, "doi_to_pmc"

    return src, "source_url_fallback"


def audit_entry(entry: dict) -> dict:
    """Returns a 'provenance' dict to be merged onto the entry."""
    gene = (entry.get("gene_name") or "").strip()
    uacc_field = (entry.get("uniprot_accession") or "").strip()
    ncbi_field = (entry.get("ncbi_or_refseq_locus") or "").strip()
    organism = (entry.get("organism") or "").strip()
    url, url_label = resolve_paper_url(entry)

    # Short-circuit the already-known GUESS cases
    if "community ortholog" in uacc_field.lower() or "community ortholog" in ncbi_field.lower():
        return {
            "verdict": "GUESS",
            "rationale": "JSON itself labels this as a community ortholog",
            "uniprot_in_paper": False,
            "ncbi_in_paper": False,
            "candidate_accessions_near_gene": [],
            "fetch_status": "skipped",
        }
    if "(family-level" in ncbi_field or "PF" in organism:
        return {
            "verdict": "GUESS",
            "rationale": "JSON specifies only a Pfam family, no locus",
            "uniprot_in_paper": False,
            "ncbi_in_paper": False,
            "candidate_accessions_near_gene": [],
            "fetch_status": "skipped",
        }

    if not url:
        return {
            "verdict": "UNVERIFIABLE",
            "rationale": "No source_url in JSON entry",
            "uniprot_in_paper": False,
            "ncbi_in_paper": False,
            "candidate_accessions_near_gene": [],
            "fetch_status": "no_url",
        }

    print(f"  fetching {url} ...")
    html, status = fetch_url(url)
    time.sleep(FETCH_DELAY)
    if not html or status != 200:
        return {
            "verdict": "UNVERIFIABLE",
            "rationale": f"Fetch failed (HTTP status={status})",
            "uniprot_in_paper": False,
            "ncbi_in_paper": False,
            "candidate_accessions_near_gene": [],
            "fetch_status": f"http_{status}",
        }

    text = strip_html(html)

    uacc_clean = extract_bare_accession(uacc_field)
    ncbi_clean = extract_bare_accession(ncbi_field)

    uacc_in = bool(uacc_clean and uacc_clean.lower() in text.lower())
    ncbi_in = bool(ncbi_clean and ncbi_clean.lower() in text.lower())

    # Look for accessions near the gene name regardless of what was already claimed
    gene_first_token = re.split(r"[ /(_]", gene)[0] if gene else ""
    candidates = near_gene_accessions(text, gene_first_token) if gene_first_token else []
    # de-duplicate against the already-claimed accession
    candidates = [c for c in candidates if c not in {uacc_clean, ncbi_clean}]

    if uacc_in or ncbi_in:
        verdict = "IN_PAPER"
        what = []
        if uacc_in:
            what.append(f"UniProt {uacc_clean}")
        if ncbi_in:
            what.append(f"NCBI {ncbi_clean}")
        rationale = "accession(s) literally appear in fetched paper text: " + ", ".join(what)
    elif uacc_clean or ncbi_clean:
        verdict = "INFERRED"
        rationale = (
            f"claimed accession(s) not found in fetched paper text; "
            f"likely looked up via UniProt/NCBI search by gene+organism. "
            f"Checked: {uacc_clean or ''} {ncbi_clean or ''}".strip()
        )
    else:
        verdict = "UNVERIFIABLE"
        rationale = "JSON has no parseable accession in either field"

    return {
        "verdict": verdict,
        "rationale": rationale,
        "uniprot_field_raw": uacc_field,
        "uniprot_clean": uacc_clean,
        "uniprot_in_paper": uacc_in,
        "ncbi_field_raw": ncbi_field,
        "ncbi_clean": ncbi_clean,
        "ncbi_in_paper": ncbi_in,
        "candidate_accessions_near_gene": candidates,
        "url_fetched": url,
        "url_label": url_label,
        "fetch_status": "ok",
    }


def main() -> int:
    entries = json.loads(SRC.read_text())
    audited: list[dict] = []

    print(f"Auditing {len(entries)} protein records...\n")
    for i, e in enumerate(entries, 1):
        gene = e.get("gene_name", "?")
        print(f"[{i}/{len(entries)}] {gene}")
        prov = audit_entry(e)
        merged = dict(e)
        merged["provenance"] = prov
        audited.append(merged)
        print(f"  -> {prov['verdict']}: {prov['rationale'][:160]}")
        if prov.get("candidate_accessions_near_gene"):
            print(f"     other accessions near gene name: "
                  f"{prov['candidate_accessions_near_gene'][:6]}")
        print()

    OUT_JSON.write_text(json.dumps(audited, indent=2))

    # TSV report
    cols = [
        "gene_name", "organism", "verdict",
        "uniprot_claimed", "uniprot_in_paper",
        "ncbi_claimed", "ncbi_in_paper",
        "candidate_accessions_near_gene", "rationale", "fetch_status",
    ]
    with OUT_TSV.open("w") as fh:
        fh.write("\t".join(cols) + "\n")
        for e in audited:
            p = e["provenance"]
            row = [
                e.get("gene_name", ""),
                e.get("organism", ""),
                p["verdict"],
                p.get("uniprot_clean", "") or "",
                "yes" if p.get("uniprot_in_paper") else "no",
                p.get("ncbi_clean", "") or "",
                "yes" if p.get("ncbi_in_paper") else "no",
                ";".join(p.get("candidate_accessions_near_gene", []) or []),
                p["rationale"],
                p.get("fetch_status", ""),
            ]
            fh.write("\t".join(str(c).replace("\t", " ") for c in row) + "\n")

    # Summary
    print("\n=== SUMMARY ===")
    from collections import Counter
    c = Counter(a["provenance"]["verdict"] for a in audited)
    for k in ("IN_PAPER", "INFERRED", "GUESS", "UNVERIFIABLE"):
        print(f"  {k:>14s}: {c.get(k, 0)}")
    print(f"\nWrote: {OUT_JSON.name}, {OUT_TSV.name}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
