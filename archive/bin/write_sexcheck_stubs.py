#!/usr/bin/env python3
"""
write_sexcheck_stubs.py
Writes empty sex-check output files when the process cannot run.
"""
import argparse


def write_tsv(path, header):
    with open(path, "w", encoding="utf-8") as handle:
        handle.write(header + "\n")


def write_sexcheck_stubs(out_dir="."):
    write_tsv(
        f"{out_dir}/sexcheck.txt",
        "FID\tIID\tPEDSEX\tSNPSEX\tSTATUS\tF\tCOLLECTED_SEX\tINFERRED_SEX\tDISCORDANT",
    )
    write_tsv(
        f"{out_dir}/sex_discordance.txt",
        "FID\tIID\tCOLLECTED_SEX\tINFERRED_SEX\tF\tSTATUS",
    )
    write_tsv(f"{out_dir}/sexcheck_multimind.tsv", "IID\tmind_base")
    write_tsv(
        f"{out_dir}/sexcheck_plate_report.tsv",
        "PLATE\tN_SAMPLES\tN_DISCORDANT\tPCT_DISCORDANT",
    )


def main():
    parser = argparse.ArgumentParser(description="Write empty sex-check stub outputs")
    parser.add_argument("--reason", default="Sex check unavailable")
    args = parser.parse_args()

    write_sexcheck_stubs()
    print(f"[write_sexcheck_stubs] {args.reason}", flush=True)


if __name__ == "__main__":
    main()
