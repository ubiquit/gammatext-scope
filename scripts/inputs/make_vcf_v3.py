"""Build the per-chromosome normal and tumour VCFs for the melanoma v3 consensus.

normal_v3  = germline records (file order preserved) minus records whose
             reference span overlaps a somatic record's span
tumour_v3  = normal_v3 with the somatic records merged in by position
             (germline records first at equal positions)

Both files therefore present the germline records in the same order, so
`bcftools consensus` makes the same apply/skip decisions in both samples and
the two consensus genomes differ only at the somatic records.

usage: make_vcf_v3.py <germline.vcf.gz> <somatic.vcf.gz> <normal_out.vcf> <tumour_out.vcf>
"""
import bisect, gzip, sys

germ_path, som_path, out_n, out_t = sys.argv[1:5]

def read(path):
    header, recs = [], []
    with gzip.open(path, 'rt', errors='ignore') as fh:
        for line in fh:
            if line.startswith('#'): header.append(line)
            else: recs.append(line.rstrip('\n'))
    return header, recs

g_header, g_recs = read(germ_path)
_, s_recs = read(som_path)

def pos_ref(rec):
    f = rec.split('\t', 5)
    return int(f[1]), f[3]

# somatic spans (1-based, inclusive), sorted
spans = sorted((pos, pos + len(ref) - 1) for pos, ref in map(pos_ref, s_recs))
starts = [s for s, _ in spans]
def overlaps_somatic(pos, ref):
    end = pos + len(ref) - 1
    i = bisect.bisect_right(starts, end)          # spans starting at or before `end`
    for s, e in spans[max(0, i - 3):i]:
        if s <= end and e >= pos: return True
    return False

kept, dropped = [], 0
for rec in g_recs:
    pos, ref = pos_ref(rec)
    if overlaps_somatic(pos, ref): dropped += 1
    else: kept.append(rec)

# somatic records re-formatted to FORMAT=GT with the record's genotype
def to_gt_only(rec):
    f = rec.split('\t')
    gt = f[9].split(':')[0] if len(f) > 9 else '0/1'
    return '\t'.join(f[:8] + ['GT', gt])
som = sorted(((pos_ref(r)[0], to_gt_only(r)) for r in s_recs), key=lambda x: x[0])

with open(out_n, 'w', newline='\n') as fn, open(out_t, 'w', newline='\n') as ft:
    for h in g_header: fn.write(h); ft.write(h)
    j = 0
    for rec in kept:
        pos = pos_ref(rec)[0]
        while j < len(som) and som[j][0] < pos:        # somatic records strictly before this germline position
            ft.write(som[j][1] + '\n'); j += 1
        fn.write(rec + '\n'); ft.write(rec + '\n')
    while j < len(som): ft.write(som[j][1] + '\n'); j += 1

print(f"germline {len(g_recs)} -> kept {len(kept)} (dropped {dropped} overlapping somatic spans); somatic {len(s_recs)}; union {len(kept) + len(s_recs)}")
