# ingest.py
# Build a local RAG index: parse -> chunk -> embed -> store/
import os, sys, json, time, glob, hashlib, argparse
from typing import List, Dict, Any
import re


# --- env & OpenAI ---
from dotenv import load_dotenv, find_dotenv
load_dotenv(find_dotenv())
from openai import OpenAI
OPENAI_API_KEY = os.getenv("OPENAI_API_KEY")
if not OPENAI_API_KEY:
    raise SystemExit("❌ OPENAI_API_KEY not found (set it in .env).")
client = OpenAI(api_key=OPENAI_API_KEY)

# --- parsing (Unstructured) ---
from unstructured.partition.pdf import partition_pdf
from unstructured.partition.docx import partition_docx

# --- tokenization for chunk sizing ---
import tiktoken
ENC = tiktoken.get_encoding("cl100k_base")

# --- FAISS vector store ---
import faiss
import numpy as np

SUPPORTED = (".pdf", ".docx")

def sha256(path: str) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for chunk in iter(lambda: f.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()

def discover_files(data_dir: str) -> List[str]:
    files = []
    for ext in SUPPORTED:
        files += [p for p in glob.glob(os.path.join(data_dir, "**", f"*{ext}"), recursive=True)]
    return sorted(files)

def parse_file(path: str, use_ocr: bool) -> List[Dict[str, Any]]:
    """Return list of blocks with text + category + page (when available)."""
    if path.lower().endswith(".pdf"):
        kwargs = {}
        if use_ocr:
            kwargs = {"strategy": "hi_res", "ocr_languages": "eng"}
        elements = partition_pdf(filename=path, **kwargs)
    elif path.lower().endswith(".docx"):
        elements = partition_docx(filename=path)
    else:
        return []

    blocks = []
    for e in elements:
        txt = (getattr(e, "text", "") or "").strip()
        if not txt:
            continue
        cat = getattr(e, "category", "Text")
        meta = getattr(e, "metadata", None)
        page = getattr(meta, "page_number", None) if meta else None
        blocks.append({"text": txt, "category": cat, "page": page})
    return blocks

def ntoks(s: str) -> int:
    return len(ENC.encode(s))

def chunk_blocks(blocks: List[Dict[str, Any]],
                 target_min=450, target_max=700, overlap=90) -> List[Dict[str, Any]]:
    """
    Section-aware chunking:
    - Titles start a new section
    - Keep chunks within a section
    - 10–15% token overlap between adjacent chunks (same section)
    """
    chunks = []
    section = None
    cur_texts, cur_pages = [], []
    cur_tokens = 0

    def flush():
        nonlocal cur_texts, cur_pages, cur_tokens, section
        if not cur_texts:
            return
        full = "\n".join(cur_texts).strip()
        if not full:
            cur_texts, cur_pages, cur_tokens = [], [], 0
            return

        # NEW: filter out None pages before taking min
        pages_num = [p for p in cur_pages if p is not None]
        page_val = min(pages_num) if pages_num else None

        chunks.append({
            "text": full,
            "section": section,
            "page": page_val,
        })
        cur_texts, cur_pages, cur_tokens = [], [], 0


    for b in blocks:
        cat = b["category"]
        txt = b["text"]

        if cat == "Title":
            # end previous section
            flush()
            section = txt
            # skip references sections completely
            if section and "reference" in section.lower():
                section = "References"
            continue

        if section and "reference" in section.lower():
            section = "References"

        t = ntoks(txt)
        # if adding would exceed target_max, flush with overlap
        if cur_texts and cur_tokens + t > target_max:
            if overlap > 0:
                joined = "\n".join(cur_texts)
                toks = ENC.encode(joined)
                tail = ENC.decode(toks[max(0, len(toks)-overlap):])
                flush()
                if tail:
                    cur_texts = [tail]
                    cur_pages = [cur_pages[-1]] if cur_pages else []
                    cur_tokens = ntoks(tail)
            else:
                flush()

        cur_texts.append(txt)
        cur_pages.append(b.get("page"))
        cur_tokens += t

    flush()
    # Merge very short chunks with their neighbors when possible
    merged = []
    for ch in chunks:
        if merged and ntoks(ch["text"]) < target_min and merged[-1]["section"] == ch["section"]:
            merged[-1]["text"] += "\n" + ch["text"]
            merged[-1]["page"] = merged[-1]["page"] or ch["page"]
        else:
            merged.append(ch)
    # Final filter: drop empty
    return [c for c in merged if c["section"] != "__DROP__" and c["text"].strip()]

def embed_texts(texts: List[str], model: str, batch_size: int = 96) -> np.ndarray:
    vecs = []
    for i in range(0, len(texts), batch_size):
        batch = texts[i:i+batch_size]
        resp = client.embeddings.create(model=model, input=batch)
        vecs.extend([d.embedding for d in resp.data])
        time.sleep(0.05)  # be polite
    arr = np.array(vecs, dtype="float32")
    # normalize for cosine via inner product
    norms = np.linalg.norm(arr, axis=1, keepdims=True) + 1e-12
    return arr / norms

def load_or_init_index(store_dir: str, dim: int) -> faiss.Index:
    path = os.path.join(store_dir, "index.faiss")
    if os.path.exists(path):
        idx = faiss.read_index(path)
        # ensure dimension matches
        if idx.d != dim:
            raise SystemExit(f"❌ Existing index dim {idx.d} != embedding dim {dim}")
        return idx
    # FlatIP because we normalized vectors (cosine)
    return faiss.IndexFlatIP(dim)

def save_index(index: faiss.Index, store_dir: str):
    os.makedirs(store_dir, exist_ok=True)
    faiss.write_index(index, os.path.join(store_dir, "index.faiss"))

def append_catalog(store_dir: str, rows: List[Dict[str, Any]]):
    os.makedirs(store_dir, exist_ok=True)
    with open(os.path.join(store_dir, "catalog.jsonl"), "a", encoding="utf-8") as f:
        for r in rows:
            f.write(json.dumps(r, ensure_ascii=False) + "\n")

def load_seen(store_dir: str) -> Dict[str, str]:
    p = os.path.join(store_dir, "seen.json")
    if os.path.exists(p):
        return json.load(open(p))
    return {}

def save_seen(store_dir: str, seen: Dict[str, str]):
    with open(os.path.join(store_dir, "seen.json"), "w") as f:
        json.dump(seen, f, indent=2)

DOI_RE = re.compile(r"\b10\.\d{4,9}/[-._;()/:A-Z0-9]+\b", re.I)
YEAR_RE = re.compile(r"\b(19|20)\d{2}\b")

def extract_citation(blocks: List[Dict[str, Any]], filename: str) -> Dict[str, Any]:
    """
    Very lightweight heuristic:
    - title: first Title block (if any), else filename sans extension
    - authors: next 1–3 short lines with commas/initials
    - year: first 19xx/20xx seen near the beginning
    - doi: first DOI-like string
    - url: DOI link if DOI found
    """
    title = None
    authors: List[str] = []
    year = None
    doi = None

    # Look at the first ~40 blocks for metadata-ish lines
    head = blocks[:40]

    # Title
    for b in head:
        if b.get("category") == "Title" and b.get("text"):
            t = b["text"].strip()
            if len(t.split()) >= 2:
                title = t
                break
    if not title:
        title = os.path.splitext(os.path.basename(filename))[0]

    # DOI and year scan
    for b in head:
        txt = b.get("text") or ""
        if not doi:
            m = DOI_RE.search(txt)
            if m:
                doi = m.group(0)
        if not year:
            y = YEAR_RE.search(txt)
            if y:
                try:
                    yr = int(y.group(0))
                    if 1900 <= yr <= 2100:
                        year = yr
                except:
                    pass
        if doi and year:
            break

    # Authors heuristic: short lines with commas/initials (after title)
    # e.g., "Doe J., Roe A., Smith B."
    after_title_started = False if title else True
    for b in head:
        txt = (b.get("text") or "").strip()
        if not txt:
            continue
        if not after_title_started and txt == title:
            after_title_started = True
            continue
        if after_title_started:
            # short-ish line with commas or semicolons and letters/dots
            if (len(txt) <= 200) and ("," in txt or ";" in txt) and any(ch.isalpha() for ch in txt):
                # split by comma/semicolon, strip
                parts = [p.strip() for p in re.split(r"[;,]", txt) if p.strip()]
                # crude filter: keep tokens that look like names (at least one space or an initial)
                candidates = [p for p in parts if (" " in p) or re.search(r"\b[A-Z]\.", p)]
                if candidates:
                    authors = candidates[:12]  # cap to avoid bloat
                    break

    url = f"https://doi.org/{doi}" if doi else None

    return {
        "title": title,
        "authors": authors or None,
        "journal": None,  # unknown unless you add a smarter parser later
        "year": year,
        "doi": doi,
        "url": url,
    }


def main():
    ap = argparse.ArgumentParser(description="Ingest PDFs/DOCX into a FAISS RAG index")
    ap.add_argument("--data", default="data", help="folder with source files")
    ap.add_argument("--store", default="rag_store", help="output folder for index/catalog")
    ap.add_argument("--model", default="text-embedding-3-large", help="OpenAI embedding model")
    ap.add_argument("--rebuild", action="store_true", help="ignore seen.json and rebuild everything")
    ap.add_argument("--ocr", action="store_true", help="enable OCR for scanned PDFs")
    args = ap.parse_args()

    files = discover_files(args.data)
    if not files:
        print(f"🤷 No PDF/DOCX found under {args.data}")
        return

    print(f"🔎 Found {len(files)} file(s). OCR={'on' if args.ocr else 'off'}")

    seen = {} if args.rebuild else load_seen(args.store)

    # Probe embedding dim
    probe = client.embeddings.create(model=args.model, input=["probe"]).data[0].embedding
    dim = len(probe)
    index = load_or_init_index(args.store, dim)

    added_points = 0
    catalog_rows = []
    ids = []

    for path in files:
        fp = sha256(path)
        if seen.get(path) == fp and not args.rebuild:
            print(f"⏩ Skipping unchanged: {os.path.basename(path)}")
            continue

        print(f"📄 Parsing: {os.path.basename(path)}")
        blocks = parse_file(path, use_ocr=args.ocr)
        citation = extract_citation(blocks, os.path.basename(path))
        if not blocks:
            print(f"   ⚠️  No text extracted; consider --ocr if PDF is scanned.")
            continue

        chunks = chunk_blocks(blocks)
        if not chunks:
            print("   ⚠️  Parser returned empty chunks.")
            continue

        texts = [c["text"] for c in chunks]
        print(f"   ✂️  {len(chunks)} chunk(s)")

        vecs = embed_texts(texts, model=args.model)
        print(f"   🔢 Embedded: {vecs.shape}")

        # add to FAISS
        index.add(vecs)
        # assign IDs as current index size minus added batch
        start_id = ids[-1] + 1 if ids else 0
        batch_ids = list(range(start_id, start_id + vecs.shape[0]))
        ids.extend(batch_ids)

        # write catalog rows
        for local_i, ch in enumerate(chunks):
            catalog_rows.append({
                "id": batch_ids[local_i],
                "paper": os.path.basename(path),
                "section": ch.get("section"),
                "page": ch.get("page"),
                "text": ch["text"],
                "source_citation": citation
            })

        seen[path] = fp
        added_points += vecs.shape[0]

        # flush periodically to disk
        if added_points >= 256:
            save_index(index, args.store)
            append_catalog(args.store, catalog_rows)
            save_seen(args.store, seen)
            catalog_rows, added_points = [], 0

    # final flush
    save_index(index, args.store)
    if catalog_rows:
        append_catalog(args.store, catalog_rows)
    save_seen(args.store, seen)

    print("✅ Done.")
    print(f"   Index: {os.path.join(args.store,'index.faiss')}")
    print(f"   Catalog: {os.path.join(args.store,'catalog.jsonl')}")
    print(f"   Seen: {os.path.join(args.store,'seen.json')}")

if __name__ == "__main__":
    main()
