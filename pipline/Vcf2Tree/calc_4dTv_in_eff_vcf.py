#!/usr/bin/env python3
# Auhtor:徐州更
from sys import argv
from pysam import VariantFile
from pysam import FastaFile

file_in = argv[1]
file_out = argv[2]
fafile = argv[3]

codon = set(["TC", "CT", "CC", "CG", "AC", "GT", "GC", "GG"])
rev_dict = dict(A='T',T='A', C='G', G='C')

bcf_in = VariantFile(file_in)
bcf_out = VariantFile(file_out, "w", header = bcf_in.header)
fa_in = FastaFile(fafile)

for rec in bcf_in.fetch():
    ann = rec.info['ANN']
    info = rec.info['ANN'][0].split('|')
    # only use synonymouse variants
    if info[1] != "synonymous_variant":
        continue
    # only the 3rd position can be 4dTv
    if int(info[9][2:-3]) % 3 != 0:
        continue
    # determine the strand by the REF column and mutation
    # if the ref is not same as the mutation site
    if rec.ref == info[9][-3]:
        pre = fa_in.fetch(rec.chrom, rec.pos-3, rec.pos-1)
    else:
        tmp = fa_in.fetch(rec.chrom, rec.pos, rec.pos+2)
        tmp.upper()
        pre = rev_dict[tmp[1]] + rev_dict[tmp[0]]
    if pre not in codon:
        continue
    bcf_out.write(rec)