#!/usr/bin/env python3
"""
rename_chr23.py
Remaps chromosome 23 to X in a PLINK BIM file so that --split-x b37 works.

Usage: python3 rename_chr23.py <input.bim> <output.bim>
"""
import sys

infile  = sys.argv[1]
outfile = sys.argv[2]

with open(infile) as fin, open(outfile, "w") as fout:
    for line in fin:
        parts = line.split("\t")
        if parts[0].strip() == "23":
            parts[0] = "X"
        fout.write("\t".join(parts))
