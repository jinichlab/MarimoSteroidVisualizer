"""Rewrite new_proteins_18.fasta with clean accession-only headers.

The Rostlab ProtT5 embedder uses the entire FASTA header line as the
record key in the output h5. Without normalization we'd get keys like
'>tr|Q5LF84|Q5LF84_BACFN Choloylglycine hydrolase OS=...' which are
unusable as join keys against protein_sequence_embedding.csv. Reduce
each header to the bare accession (Q5LF84, MFU7516964_1, etc.).
"""
import re
from pathlib import Path

HERE = Path(__file__).resolve().parent
SRC = HERE / "new_proteins_18.fasta"
OUT = HERE / "new_proteins_18.clean.fasta"

ACC_FROM_TR = re.compile(r"^(?:sp|tr)\|([A-Z0-9]+)\|")

def normalize_header(raw: str) -> str:
    s = raw.lstrip(">").strip()
    m = ACC_FROM_TR.match(s)
    if m:
        return m.group(1)
    # NCBI-style "MFU7516964.1 ..." — take token before space, replace . with _
    token = s.split()[0]
    return token.replace(".", "_").replace("/", "_")

def main() -> None:
    out_lines: list[str] = []
    with SRC.open() as fh:
        for line in fh:
            if line.startswith(">"):
                acc = normalize_header(line)
                out_lines.append(f">{acc}\n")
            else:
                out_lines.append(line)
    OUT.write_text("".join(out_lines))
    print(f"Wrote {OUT}")
    print("Accessions:")
    for ln in out_lines:
        if ln.startswith(">"):
            print(f"  {ln.strip()}")

if __name__ == "__main__":
    main()
