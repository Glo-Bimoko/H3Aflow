#!/usr/bin/env python3
"""
make_sex_update.py
Converts sex_info.tsv (sampleid, sex: 0=Female 1=Male) to a PLINK
--update-sex file (FID  IID  SEX: 1=Male 2=Female 0=Unknown).

Usage: python3 make_sex_update.py <sex_info.tsv> <plink_sex_update.txt>
"""
import sys

infile  = sys.argv[1]
outfile = sys.argv[2]

with open(infile) as fh, open(outfile, "w") as out:
    fh.readline()  # skip header
    for line in fh:
        parts = line.strip().split("\t")
        if len(parts) < 2:
            continue
        iid  = parts[0].strip()
        code = parts[1].strip()
        plink_sex = "1" if code == "1" else ("2" if code == "0" else "0")
        out.write("0\t" + iid + "\t" + plink_sex + "\n")
