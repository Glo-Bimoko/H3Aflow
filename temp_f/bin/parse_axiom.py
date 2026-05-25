#!/usr/bin/env python3
import sys
import pandas as pd
import re
import tempfile
from pathlib import Path

def main():
    if len(sys.argv) != 3:
        print("Usage: parse_axiom.py <input_file> <output_file>")
        sys.exit(1)

    infile = Path(sys.argv[1])
    outfile = Path(sys.argv[2])

    GENO_COL_RE = re.compile(r"^(.+?)\.CEL_call_code$", re.IGNORECASE)

    try:
        print(f"Processing file: {infile}")

        if infile.suffix.lower() in [".xlsx", ".xls"]:
            raw = pd.read_excel(infile, header=None, dtype=str)
        else:
            with open(infile, 'r', encoding='utf-8') as f:
                lines = f.readlines()

            print(f"File has {len(lines)} lines")

            data_start = 0
            for i, line in enumerate(lines):
                if not line.strip().startswith('##'):
                    data_start = i
                    print(f"Data starts at line {i}: {line.strip()[:100]}...")
                    break

            if data_start >= len(lines) - 1:
                print("ERROR: No data found after metadata lines")
                sys.exit(1)

            with tempfile.NamedTemporaryFile(mode='w', suffix='.txt', delete=False) as tmp_file:
                for line in lines[data_start:]:
                    tmp_file.write(line)
                tmp_filename = tmp_file.name

            try:
                raw = pd.read_csv(tmp_filename, sep='\t', dtype=str, on_bad_lines='warn')
            except TypeError:
                raw = pd.read_csv(tmp_filename, sep='\t', dtype=str,
                                  error_bad_lines=False, warn_bad_lines=False)
            except Exception as e:
                print(f"Error reading with pandas: {e}")
                with open(tmp_filename, 'r') as f:
                    clean_lines = []
                    header = None
                    for line_num, line in enumerate(f):
                        line = line.strip()
                        if not line:
                            continue
                        parts = line.split('\t')
                        if header is None:
                            header = parts
                            expected_cols = len(parts)
                        elif len(parts) == expected_cols:
                            clean_lines.append(parts)
                        else:
                            print(f"Skipping malformed line {line_num + 1}: "
                                  f"expected {expected_cols} fields, got {len(parts)}")
                if header and clean_lines:
                    raw = pd.DataFrame(clean_lines, columns=header)
                else:
                    raise Exception("No valid data found")
            finally:
                Path(tmp_filename).unlink()

        print(f"Data shape: {raw.shape}")
        print(f"Columns: {list(raw.columns)}")

        # --- Detect array type from annotation column ---
        array_type = "unknown"
        # Re-read metadata lines to find annotation file
        with open(infile, 'r', encoding='utf-8') as f:
            for line in f:
                if line.startswith('##annotation-file'):
                    if 'PMDA' in line:
                        array_type = 'PMDA'
                    elif '3x4' in line or '3X4' in line:
                        array_type = '3x4'
                    break
        print(f"Array type detected: {array_type}")

        # --- Find rsID column ---
        # Try exact match first, then case-insensitive
        rsid_col = None
        for candidate in ["extended_rsid", "rsid", "RS ID", "rsID", "snp_id"]:
            if candidate in raw.columns:
                rsid_col = candidate
                break
        if rsid_col is None:
            # Case-insensitive fallback
            for col in raw.columns:
                if 'rsid' in col.lower() or 'rs_id' in col.lower():
                    rsid_col = col
                    break
        if rsid_col is None:
            print("ERROR: Could not find an rsID column.")
            print(f"Available columns: {list(raw.columns)}")
            sys.exit(1)
        print(f"Found rsID column: '{rsid_col}'")

        # --- Find genotype call column ---
        geno_cols = [c for c in raw.columns if GENO_COL_RE.match(str(c))]
        if not geno_cols:
            print("ERROR: Could not find genotype column matching '*.CEL_call_code'")
            for col in raw.columns:
                print(f"  - {col}")
            sys.exit(1)

        geno_col = geno_cols[0]
        print(f"Found genotype column: '{geno_col}'")

        match = GENO_COL_RE.match(geno_col)
        sample_id = match.group(1) if match else infile.stem.replace('.CEL', '')
        print(f"Sample ID: {sample_id}")

        # --- Build output: rsid + genotype, index on rsid ---
        result_df = raw[[rsid_col, geno_col]].copy()
        result_df = result_df.rename(columns={rsid_col: "rsid", geno_col: sample_id})

        initial_rows = len(result_df)

        # Drop rows with no rsID — these are array QC probes (e.g. AFFX-QC-*)
        # that have no biological meaning and no cross-array identity
        result_df = result_df[result_df["rsid"].notna()]
        result_df = result_df[~result_df["rsid"].str.strip().isin(["", "---", "NA", "N/A"])]
        no_rsid_dropped = initial_rows - len(result_df)
        print(f"Dropped {no_rsid_dropped} probes with no rsID (QC/control probes)")

        # Keep NaN genotype calls — don't drop them here.
        # The comparison step will skip NaN calls per pair.
        result_df = result_df.set_index("rsid")

        print(f"Final rows (SNPs with rsID): {len(result_df)}")
        result_df.to_csv(outfile, sep="\t", index=True)
        print(f"Output saved to: {outfile}")
        print(f"First few rows:\n{result_df.head()}")

    except Exception as e:
        print(f"ERROR processing {infile}: {str(e)}")
        import traceback
        traceback.print_exc()
        sys.exit(1)

if __name__ == "__main__":
    main()
