"""Triage every paper in dorrestein_catalog_2024_2026.json for the
presence of extractable novel proteins.

For each paper:
  1. Resolve DOI -> PMC ID via NCBI idconv. If no PMC, fall back to the
     PubMed abstract via E-utilities efetch (db=pubmed, rettype=abstract).
  2. Fetch the best-available text (PMC full text > PubMed abstract > nothing).
  3. Run two scans:
        - accession_patterns: UniProt-shaped, NCBI-protein-shaped, locus-tag-shaped
        - enzyme_keywords: hydrolase / hydroxylase / reductase / etc.
  4. Verdict:
        HAS_ACCESSIONS    one or more concrete accessions in the paper text
        ENZYME_LIKELY     enzyme keywords present (>=2 distinct), no accessions
        WEAK_SIGNAL       1 keyword hit
        NO_SIGNAL         no enzyme words at all
        UNVERIFIABLE      no full text + no abstract available

Writes triage_catalog_report.tsv ordered worst-known-first so you can
review the high-signal rows first.
"""
from __future__ import annotations

import json
import re
import time
import urllib.parse
import urllib.request
from collections import Counter
from pathlib import Path

HERE = Path(__file__).resolve().parent
SRC = HERE / "dorrestein_catalog_2024_2026.json"
OUT = HERE / "triage_catalog_report.tsv"
OUT_JSON = HERE / "triage_catalog_report.json"

UA = "literature-triage (adsiordia@ucsd.edu)"
FETCH_DELAY = 0.4  # NCBI E-utilities allows 3 req/sec without a key; play safe

UNIPROT_RX = re.compile(
    r"\b("
    r"[OPQ][0-9][A-Z0-9]{3}[0-9]"
    r"|[A-NR-Z][0-9](?:[A-Z][A-Z0-9]{2}[0-9]){1,2}"
    r")\b"
)
NCBI_PROT_RX = re.compile(
    r"\b("
    r"(?:WP_|NP_|YP_|XP_|MFU|MFV)\d+(?:\.\d+)?"
    r"|[A-Z]{3}\d{5,}\.\d+"
    r")\b"
)
LOCUS_RX = re.compile(
    r"\b("
    r"Elen_\d{3,5}"
    r"|BF9343_\d{3,5}"
    r"|HMPREF\d+_\d+"
    r"|ci\d{3,5}"
    r"|bai[A-Z]{1,2}"
    r"|bsh[/_]?\w?"
    r")\b"
)
ENZYME_RX = re.compile(
    r"\b("
    r"hydrolase|hydroxylase|reductase|dehydrogenase|transferase|oxidoreductase"
    r"|cytochrome P450|CYP\d+\w*|HSD\d+\w*|monooxygenase|isomerase|sulfat\w+"
    r"|methyltransferase|acyltransferase|kinase|lyase|epimerase|deconjugat\w+"
    r"|sterol\s+\w+ase|steroid\s+\w+ase|bile salt hydrolase|BSH/T|baiCD|baiE|baiH|baiI|baiA"
    r"|oxysterol-\w+|amine-conjug\w+"
    r")\b",
    re.I,
)


def http_get(url: str, timeout: int = 30) -> str:
    req = urllib.request.Request(url, headers={"User-Agent": UA})
    try:
        with urllib.request.urlopen(req, timeout=timeout) as r:
            return r.read().decode("utf-8", errors="replace")
    except Exception:
        return ""


def strip_html(html: str) -> str:
    html = re.sub(r"<script.*?</script>", " ", html, flags=re.S | re.I)
    html = re.sub(r"<style.*?</style>", " ", html, flags=re.S | re.I)
    html = re.sub(r"<[^>]+>", " ", html)
    return re.sub(r"\s+", " ", html)


def doi_to_pmc(doi: str) -> str | None:
    if not doi:
        return None
    api = (
        "https://www.ncbi.nlm.nih.gov/pmc/utils/idconv/v1.0/"
        f"?ids={urllib.parse.quote(doi)}&format=json&tool=triage&email=adsiordia@ucsd.edu"
    )
    txt = http_get(api)
    if not txt:
        return None
    try:
        for rec in json.loads(txt).get("records", []):
            if rec.get("pmcid"):
                return rec["pmcid"]
    except Exception:
        pass
    return None


def fetch_pubmed_abstract(pmid: str) -> str:
    if not pmid:
        return ""
    url = (
        "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi"
        f"?db=pubmed&id={pmid}&rettype=abstract&retmode=text"
    )
    return http_get(url)


def triage_paper(paper: dict, idx: int, total: int) -> dict:
    title = paper.get("title", "")
    doi = paper.get("doi", "") or ""
    pmid = paper.get("pmid", "") or ""

    text_source = ""
    text = ""

    if doi and doi != "unknown":
        pmcid = doi_to_pmc(doi)
        time.sleep(FETCH_DELAY)
        if pmcid:
            html = http_get(f"https://pmc.ncbi.nlm.nih.gov/articles/{pmcid}/")
            time.sleep(FETCH_DELAY)
            if html:
                text = strip_html(html)
                text_source = f"PMC:{pmcid}"
    if not text and pmid:
        abs_text = fetch_pubmed_abstract(pmid)
        time.sleep(FETCH_DELAY)
        if abs_text:
            text = abs_text
            text_source = f"PubMed_abstract:{pmid}"

    # Even if we couldn't fetch, the catalog title + relevance_note are something
    catalog_blob = " | ".join(
        str(paper.get(k, ""))
        for k in ("title", "relevance_note", "steroid_or_bile_acid_relevance")
    )

    scan_text = (text or "") + " " + catalog_blob
    uniprots = sorted(set(UNIPROT_RX.findall(scan_text)))
    ncbis = sorted(set(NCBI_PROT_RX.findall(scan_text)))
    locuses = sorted(set(LOCUS_RX.findall(scan_text)))
    enzyme_hits = Counter(m.lower() for m in ENZYME_RX.findall(scan_text))

    # Drop stop-tokens that look like accessions but aren't
    DROP = {"PMC", "DOI", "ATP", "NADP", "NADPH", "CoA"}
    uniprots = [u for u in uniprots if u not in DROP]

    # Verdict
    if uniprots or ncbis or locuses:
        verdict = "HAS_ACCESSIONS"
    elif len(enzyme_hits) >= 2:
        verdict = "ENZYME_LIKELY"
    elif len(enzyme_hits) == 1:
        verdict = "WEAK_SIGNAL"
    elif not text:
        verdict = "UNVERIFIABLE"
    else:
        verdict = "NO_SIGNAL"

    return {
        "idx": idx + 1,
        "year": paper.get("year", ""),
        "month": paper.get("month", ""),
        "journal": paper.get("journal_or_preprint_server", ""),
        "doi": doi,
        "pmid": pmid,
        "title": title,
        "catalog_flag": paper.get("introduces_novel_protein", "?"),
        "text_source": text_source or "(none)",
        "verdict": verdict,
        "uniprot_candidates": ";".join(uniprots[:8]),
        "ncbi_candidates": ";".join(ncbis[:8]),
        "locus_candidates": ";".join(locuses[:8]),
        "enzyme_keyword_hits": ";".join(f"{k}({v})" for k, v in enzyme_hits.most_common(6)),
    }


def main() -> None:
    catalog = json.load(SRC.open())["papers"]
    print(f"Triaging {len(catalog)} papers from catalog...\n")
    results = []
    for i, p in enumerate(catalog):
        title_short = (p.get("title") or "")[:80]
        print(f"[{i+1:>3}/{len(catalog)}] {title_short}")
        r = triage_paper(p, i, len(catalog))
        results.append(r)
        print(f"     -> {r['verdict']:<14s}  via {r['text_source']}  "
              f"acc={r['uniprot_candidates'] or r['ncbi_candidates'] or r['locus_candidates'] or '-'}")
    # Sort: HAS_ACCESSIONS > ENZYME_LIKELY > WEAK_SIGNAL > NO_SIGNAL > UNVERIFIABLE
    order = {"HAS_ACCESSIONS": 0, "ENZYME_LIKELY": 1, "WEAK_SIGNAL": 2, "NO_SIGNAL": 3, "UNVERIFIABLE": 4}
    results.sort(key=lambda r: (order.get(r["verdict"], 9), -int(r["year"]) if str(r["year"]).isdigit() else 0))

    # Write TSV
    cols = ["idx", "year", "month", "journal", "doi", "pmid",
            "verdict", "text_source", "catalog_flag",
            "uniprot_candidates", "ncbi_candidates", "locus_candidates",
            "enzyme_keyword_hits", "title"]
    with OUT.open("w") as fh:
        fh.write("\t".join(cols) + "\n")
        for r in results:
            row = [str(r.get(c, "")).replace("\t", " ").replace("\n", " ") for c in cols]
            fh.write("\t".join(row) + "\n")
    OUT_JSON.write_text(json.dumps(results, indent=2))

    print("\n=== SUMMARY ===")
    c = Counter(r["verdict"] for r in results)
    for v in ("HAS_ACCESSIONS", "ENZYME_LIKELY", "WEAK_SIGNAL", "NO_SIGNAL", "UNVERIFIABLE"):
        print(f"  {v:<15s}: {c.get(v, 0)}")
    print(f"\nWrote {OUT.name} (sorted by signal strength) and {OUT_JSON.name}")


if __name__ == "__main__":
    main()
