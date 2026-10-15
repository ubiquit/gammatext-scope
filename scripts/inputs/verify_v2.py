"""Verification of the v2 inputs: symmetry of allele application (melanoma) and
identity of the breast-cancer v2 sequences with v1."""
import glob, gzip, hashlib, random, sys
random.seed(7)
STREAM = "D:/Files/Stream"; CS = "D:/Cancer sequences"

def chr21(sample):
    f = glob.glob(f"{STREAM}/{sample}/chr21.*.txt")[0]
    return ''.join(l.strip() for l in open(f)).upper()

def ref21(fa, fai):
    d = {l.split('\t')[0]: l.split('\t') for l in open(fai)}
    n = '21' if '21' in d else 'chr21'
    L, off, lb, lw = int(d[n][1]), int(d[n][2]), int(d[n][3]), int(d[n][4])
    with open(fa, 'rb') as fh:
        fh.seek(off); raw = fh.read((L // lb + 1) * lw)
    return raw.replace(b'\n', b'').replace(b'\r', b'')[:L].decode().upper()

def snvs(vcf, chrom):
    het, hom = [], []
    with gzip.open(vcf, 'rt', errors='ignore') as fh:
        for line in fh:
            if line[0] == '#': continue
            p = line.split('\t')
            if p[0] != chrom: continue
            alts = p[4].split(',')
            if len(p[3]) != 1 or len(alts[0]) != 1 or len(alts) > 1: continue
            gt = p[9].split(':')[0].replace('|', '/') if len(p) > 9 else '1/1'
            (het if gt in ('0/1', '1/0') else hom if gt == '1/1' else []).append((int(p[1]), p[3], alts[0]))
    return het, hom

def base_at(cons, ref, pos, k=32, radius=4000):
    flank = ref[pos-1-k:pos-1]
    if 'N' in flank: return None
    lo = max(0, pos-1-k-radius); seg = cons[lo:pos-1-k+radius]; i = seg.find(flank)
    if i < 0 or seg.find(flank, i+1) >= 0: return None
    j = lo + i + k
    return cons[j] if j < len(cons) else None

def tally(sites, ref, cn, ct, nmax=400):
    out = {'normal': {'REF': 0, 'ALT': 0, 'other': 0}, 'tumour': {'REF': 0, 'ALT': 0, 'other': 0}}; used = 0
    for pos, r, a in random.sample(sites, min(nmax, len(sites))):
        if ref[pos-1] != r: continue
        bn, bt = base_at(cn, ref, pos), base_at(ct, ref, pos)
        if bn is None or bt is None: continue
        used += 1
        for lab, b in (('normal', bn), ('tumour', bt)):
            out[lab]['REF' if b == r else 'ALT' if b == a else 'other'] += 1
    return used, out

ok = True
suf = sys.argv[1] if len(sys.argv) > 1 else "v2"
print(f"== melanoma {suf}, chr21: allele application ==")
ref = ref21(f"{CS}/melanoma/GRCh38_chr.fixed.fa", f"{CS}/melanoma/GRCh38_chr.fixed.fa.fai")
cn, ct = chr21(f"HomoSapiens_Melanoma_Normal_{suf}"), chr21(f"HomoSapiens_Melanoma_Tumor_{suf}")
het, hom = snvs(f"{CS}/melanoma/tmp_norm_chr/21.vcf.gz", '21')
for label, sites, expect in (("heterozygous germline", het, 'ALT'), ("homozygous germline", hom, 'ALT')):
    used, t = tally(sites, ref, cn, ct)
    print(f"  {label:22s} ({used} sites): normal {t['normal']}  tumour {t['tumour']}")
    for side in ('normal', 'tumour'):
        if used and t[side][expect] / used < 0.97: ok = False; print(f"    !! {side} should carry {expect}")
shet, shom = snvs(f"{CS}/melanoma/tmp_som_chr/21.vcf.gz", '21')
used, t = tally(shet + shom, ref, cn, ct)
print(f"  {'somatic':22s} ({used} sites): normal {t['normal']}  tumour {t['tumour']}")
if used and (t['normal']['REF'] / used < 0.97 or t['tumour']['ALT'] / used < 0.97): ok = False; print("    !! somatic sites should be REF in normal, ALT in tumour")

import os
if os.path.isdir(f"{STREAM}/HomoSapiens_BreastCancer_Normal_{suf}"):
    print(f"== breast cancer: {suf} sequences vs v1 (chr21 md5) ==")
    for v1, v2 in (("HomoSapiens_BreastCancer_Normal", f"HomoSapiens_BreastCancer_Normal_{suf}"), ("HomoSapiens_BreastCancer_Tumor", f"HomoSapiens_BreastCancer_Tumor_{suf}")):
        a = hashlib.md5(chr21(v1).encode()).hexdigest(); b = hashlib.md5(chr21(v2).encode()).hexdigest()
        print(f"  {v2}: {'identical to v1' if a == b else 'DIFFERENT from v1 (unexpected)'}")
        if a != b: ok = False

print("VERIFY:", "PASS" if ok else "FAIL")
sys.exit(0 if ok else 1)
