"""Make the protein-coding gene mask of a tumour/normal pair identical.

The annotation is lifted onto each consensus genome separately, and paftools may
reject a gene in one sample only (an indel inside the gene breaks the length check).
Such a gene would be masked in one filtered sequence and not in the other, and from
there on every window of that chromosome would compare different DNA. Protein-coding
`gene` lines whose ID is present in only one of the two files are therefore removed
from the file that has them; everything else is left untouched.

usage: python inputs/pair_masks.py <normal genomic.gff> <tumour genomic.gff>
"""
import re, sys

def pc_gene_id(line):
    f = line.split("\t")
    if len(f) > 8 and f[2] == "gene" and "gene_type=protein_coding" in f[8]:
        m = re.search(r"ID=([^;]+)", f[8])
        return m.group(1) if m else None
    return None

def ids(path):
    with open(path, encoding="utf-8", errors="ignore") as fh:
        return {g for g in (pc_gene_id(l) for l in fh) if g}

a, b = sys.argv[1], sys.argv[2]
ia, ib = ids(a), ids(b)
for path, only in ((a, ia - ib), (b, ib - ia)):
    if not only:
        print(f"{path}: all {len(ids(path))} protein-coding genes are in the partner too"); continue
    with open(path, encoding="utf-8", errors="ignore") as fh: lines = fh.readlines()
    kept = [l for l in lines if pc_gene_id(l) not in only]
    with open(path, "w", encoding="utf-8", newline="") as fh: fh.writelines(kept)
    print(f"{path}: removed {len(lines) - len(kept)} protein-coding gene line(s) absent in the partner: {', '.join(sorted(only))}")
print(f"shared protein-coding genes: {len(ia & ib)}")
