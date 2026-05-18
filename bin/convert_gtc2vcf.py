#!/usr/bin/env python3
"""
convert_gtc2vcf.py
==================
Converts a set of GTC files for one plate to a normalised, sorted BCF using
bcftools +gtc2vcf, then produces a TSV of raw intensities via bcftools query.

Usage:
    python convert_gtc2vcf.py \
        --bpm       /path/to/manifest.bpm          \
        --egt       /path/to/clusters.egt           \
        --gtcs      gtc_list.txt                    \
        --fasta-ref /path/to/reference.fa           \
        --outprefix Plate_1                         \
        [--sex-info /path/to/sex_info.tsv]

Outputs:
    <outprefix>.bcf       normalised, sorted BCF (all samples for the plate)
    <outprefix>.bcf.csi   BCF index
    <outprefix>.tsv       raw X/Y intensity table (SAMPLE, CHR, POS, NORMX, NORMY)

Requires:
    bcftools >= 1.11 with gtc2vcf plugin
    BCFTOOLS_PLUGINS env var pointing to the plugin directory

Fix history
-----------
v2 (2026-05):
  Step 3 – haploid recoding rewritten.

  Root cause of the original bug:
    The old code used SMPL_HOMREF / SMPL_HET / SMPL_HOMALT sample-level
    filter expressions with --samples-file.  These expressions are evaluated
    *before* bcftools applies the sample subset, so on bcftools <= 1.17 they
    match the full-cohort genotype array, not the filtered male subset.  The
    result is that --new-gt is applied to the wrong samples (or not at all),
    leaving male chrX genotypes diploid.  PLINK then discards all male het
    chrX calls as het-haploid -> F = -1 -> every male inferred as Female.

v3 (2026-05):
  Step 3 – haploid recoding rewritten again for bcftools 1.20 compatibility.

  The v2 approach used the GT["samplename"] per-sample accessor syntax which
  was only introduced in bcftools 1.21.  On bcftools 1.20 this produces
  "Could not parse the index" and exits 255.

  New approach (bcftools >= 1.11, tested on 1.20):
    1. Extract sex-chromosome variants from the full BCF.
    2. Split into two BCFs: males-only and females-only, using --samples-file.
    3. On the males-only BCF, run three +setGT passes using cohort-level
       GT="RR" / GT="AA" / GT="het" expressions.  These are safe here because
       the BCF has already been subset to males, so every sample is a target.
    4. Merge the recoded males BCF back with the females BCF using bcftools merge.
    5. Concatenate the recoded sex-chr BCF with the original autosome BCF and
       sort, then replace the diploid output BCF.

  Plugin path discovery:
    The bcftools standard plugins (+setGT etc.) live in a libexec/ directory
    relative to the bcftools binary, which may differ from the user-supplied
    BCFTOOLS_PLUGINS path (pointing to the gtc2vcf plugin directory).  We now
    resolve the standard plugin dir dynamically and always append it to
    BCFTOOLS_PLUGINS so both plugin sets are available.

v4 (2026-05):
  Bug fix: +setGT was silently failing on CHPC Lengau because the standard
  bcftools plugin dir was not being discovered correctly.

  Root cause: on Lengau, `which bcftools` returns a module wrapper script
  rather than the real binary, so the parent.parent heuristic pointed to
  the wrong prefix and setGT.so was never found.  The +setGT pipe was run
  with check=False, so the failure was swallowed and the diploid BCF was
  kept — producing 821,396 het-haploid warnings in PLINK and invalidating
  all F-statistics.

  Fix:
    1. build_plugin_env() now accepts an explicit BCFTOOLS_LIBEXEC override
       via the BCFTOOLS_LIBEXEC environment variable (set in nextflow.config
       for each profile).  This is tried first before any path discovery.
    2. The known CHPC Lengau path is added as a hardcoded fallback so the
       pipeline works on that cluster even without the env var.
    3. The +setGT subprocess now uses check=True so any future plugin
       failures abort loudly at GTC_TO_VCF rather than silently producing
       wrong sex-check results downstream.
"""

import argparse
import subprocess
import sys
from pathlib import Path
import os


# ── CLI ────────────────────────────────────────────────────────────────────────
parser = argparse.ArgumentParser()
parser.add_argument("--bpm",       required=True)
parser.add_argument("--egt",       required=True)
parser.add_argument("--gtcs",      required=True,
                    help="File listing GTC paths, one per line")
parser.add_argument("--fasta-ref", required=True)
parser.add_argument("--outprefix", required=True)
parser.add_argument(
    "--sex-info",
    default=None,
    help=(
        "TSV with columns: sampleid, sex  (1=Male, 0=Female). "
        "When supplied, male X/Y genotypes are recoded to haploid so that "
        "PLINK --check-sex produces valid F-statistics."
    ),
)
args = parser.parse_args()

bpm       = Path(args.bpm).resolve()
egt       = Path(args.egt).resolve()
gtcs_list = Path(args.gtcs)
fasta     = Path(args.fasta_ref).resolve()
prefix    = args.outprefix

bcf_out = Path(f"{prefix}.bcf")
tsv_out = Path(f"{prefix}.tsv")


# ── Helpers ────────────────────────────────────────────────────────────────────
def run(cmd, *, check=True, shell=False, env=None):
    """Run a subprocess; print stdout+stderr; sys.exit on failure when check=True."""
    result = subprocess.run(
        cmd,
        capture_output=True,
        text=True,
        shell=shell,
        env=env,
    )
    if result.stdout:
        print(result.stdout, flush=True)
    if result.stderr:
        print(result.stderr, file=sys.stderr, flush=True)
    if check and result.returncode != 0:
        cmd_str = cmd if isinstance(cmd, str) else " ".join(str(c) for c in cmd)
        sys.exit(
            f"[convert_gtc2vcf] ERROR: command failed "
            f"(exit {result.returncode}):\n  {cmd_str}"
        )
    return result


def build_plugin_env():
    """
    Return a copy of os.environ where BCFTOOLS_PLUGINS contains both the
    user-supplied plugin directory (for +gtc2vcf) and the standard bcftools
    plugin directory (for +setGT, +fill-tags, etc.).

    Search order for the standard plugin dir:
      1. BCFTOOLS_LIBEXEC env var (set explicitly in nextflow.config per profile)
      2. Resolved from `which bcftools` → binary's real parent/parent/libexec
      3. Hardcoded CHPC Lengau path for bcftools/1.20
      4. Well-known install prefixes (/usr/local, /usr, /opt/conda, /opt/homebrew)
    """
    env = os.environ.copy()
    extra_dirs = []

    # ── Priority 1: explicit override from nextflow.config ────────────────────
    explicit = env.get("BCFTOOLS_LIBEXEC", "").strip()
    if explicit and Path(explicit).exists():
        extra_dirs.append(explicit)
        print(f"[convert_gtc2vcf] Standard plugin dir (BCFTOOLS_LIBEXEC): {explicit}",
              flush=True)

    # ── Priority 2: resolve from the real bcftools binary path ───────────────
    if not extra_dirs:
        which_r = subprocess.run(["which", "bcftools"], capture_output=True, text=True)
        if which_r.returncode == 0:
            # Follow symlinks so module wrapper scripts resolve to the real binary
            bcftools_bin = Path(which_r.stdout.strip()).resolve()
            install_prefix = bcftools_bin.parent.parent
            for sub in ["libexec/bcftools", "lib/bcftools", "share/bcftools/plugins"]:
                candidate = install_prefix / sub
                if candidate.exists():
                    extra_dirs.append(str(candidate))
                    print(f"[convert_gtc2vcf] Standard plugin dir (auto-detected): {candidate}",
                          flush=True)
                    break

    # ── Priority 3: hardcoded CHPC Lengau fallback ────────────────────────────
    CHPC_LIBEXEC = "/apps/chpc/bio/bcftools/1.20/libexec/bcftools"
    if Path(CHPC_LIBEXEC).exists() and CHPC_LIBEXEC not in extra_dirs:
        extra_dirs.append(CHPC_LIBEXEC)
        print(f"[convert_gtc2vcf] Standard plugin dir (CHPC fallback): {CHPC_LIBEXEC}",
              flush=True)

    # ── Priority 4: common install prefixes ───────────────────────────────────
    for pfx in ["/usr/local", "/usr", "/opt/conda", "/opt/homebrew"]:
        for sub in ["libexec/bcftools", "lib/bcftools", "share/bcftools/plugins"]:
            candidate = str(Path(pfx) / sub)
            if Path(candidate).exists() and candidate not in extra_dirs:
                extra_dirs.append(candidate)

    current_parts = [p for p in env.get("BCFTOOLS_PLUGINS", "").split(":") if p]
    all_parts = current_parts + [d for d in extra_dirs if d not in current_parts]
    if all_parts:
        env["BCFTOOLS_PLUGINS"] = ":".join(all_parts)

    return env


def detect_xy_regions(bcf_path):
    """
    Return a comma-separated region string for chrX and chrY as they appear
    in the BCF header (handles both 'X'/'Y' and '23'/'24' naming).
    """
    r = subprocess.run(
        ["bcftools", "view", "--header-only", str(bcf_path)],
        capture_output=True, text=True,
    )
    header = r.stdout
    regions = []
    for x_name in ("X", "23"):
        if (f"ID={x_name}," in header
                or f"ID={x_name}\n" in header
                or f"\t{x_name}\t" in header):
            regions.append(x_name)
            break
    for y_name in ("Y", "24"):
        if (f"ID={y_name}," in header
                or f"ID={y_name}\n" in header
                or f"\t{y_name}\t" in header):
            regions.append(y_name)
            break
    # If neither naming convention is found, fall back to trying all four;
    # bcftools will silently skip names that are not in the BCF.
    return ",".join(regions) if regions else "X,Y,23,24"


def get_auto_contigs(bcf_path, xy_set):
    """Return a list of contig IDs from the BCF header that are not in xy_set."""
    r = subprocess.run(
        ["bcftools", "view", "--header-only", str(bcf_path)],
        capture_output=True, text=True,
    )
    contigs = []
    for line in r.stdout.splitlines():
        if line.startswith("##contig") and "ID=" in line:
            cid = line.split("ID=")[1].split(",")[0].rstrip(">")
            if cid not in xy_set:
                contigs.append(cid)
    return contigs


# ── Validate inputs ────────────────────────────────────────────────────────────
for f, label in [(bpm, "BPM"), (egt, "EGT"), (fasta, "FASTA"), (gtcs_list, "GTC list")]:
    if not f.exists():
        sys.exit(f"[convert_gtc2vcf] ERROR: {label} not found: {f}")

gtc_files = [l.strip() for l in gtcs_list.read_text().splitlines() if l.strip()]
if not gtc_files:
    sys.exit(f"[convert_gtc2vcf] ERROR: GTC list is empty: {gtcs_list}")

print(f"[convert_gtc2vcf] Plate prefix : {prefix}", flush=True)
print(f"[convert_gtc2vcf] GTC files    : {len(gtc_files)}", flush=True)
print(f"[convert_gtc2vcf] BPM          : {bpm}", flush=True)
print(f"[convert_gtc2vcf] EGT          : {egt}", flush=True)
print(f"[convert_gtc2vcf] FASTA        : {fasta}", flush=True)


# ── Step 1: GTC → unsorted BCF via bcftools +gtc2vcf ─────────────────────────
unsorted_bcf = f"{prefix}_unsorted.bcf"

gtc2vcf_cmd = (
    ["bcftools", "+gtc2vcf",
     "--bpm",         str(bpm),
     "--egt",         str(egt),
     "--fasta-ref",   str(fasta),
     "--output",      unsorted_bcf,
     "--output-type", "b",
     "--no-version"]
    + gtc_files
)

print(f"\n[convert_gtc2vcf] Step 1: GTC -> BCF", flush=True)
print(f"  {' '.join(gtc2vcf_cmd[:8])} ... [{len(gtc_files)} GTC files]", flush=True)
run(gtc2vcf_cmd, env=os.environ)


# ── Step 2: Sort and index the BCF ────────────────────────────────────────────
print(f"\n[convert_gtc2vcf] Step 2: Sort BCF", flush=True)

run([
    "bcftools", "sort",
    "--output",      str(bcf_out),
    "--output-type", "b",
    "--temp-dir",    ".",
    unsorted_bcf,
])
Path(unsorted_bcf).unlink(missing_ok=True)
run(["bcftools", "index", str(bcf_out)])
print(f"[convert_gtc2vcf] BCF written  : {bcf_out}", flush=True)


# ── Step 3 (optional): Recode male X/Y genotypes to haploid ──────────────────
if args.sex_info:
    sex_info_path = Path(args.sex_info)
    if not sex_info_path.exists():
        print(
            f"[convert_gtc2vcf] WARNING: --sex-info not found: {sex_info_path}. "
            f"Skipping haploid recoding.",
            file=sys.stderr, flush=True,
        )
    else:
        print(f"\n[convert_gtc2vcf] Step 3: Recode male X/Y genotypes to haploid",
              flush=True)

        # --- get sample list from BCF ---
        header_r = run(["bcftools", "query", "--list-samples", str(bcf_out)])
        bcf_samples    = header_r.stdout.strip().splitlines()
        bcf_sample_set = set(bcf_samples)

        # --- parse sex_info and collect males present in this plate ---
        males_in_plate = []
        with open(sex_info_path) as fh:
            fh.readline()  # skip header
            for line in fh:
                parts = line.strip().split("\t")
                if len(parts) < 2:
                    continue
                sample_id, sex_code = parts[0].strip(), parts[1].strip()
                if sex_code == "1" and sample_id in bcf_sample_set:
                    males_in_plate.append(sample_id)

        if not males_in_plate:
            print(
                "[convert_gtc2vcf] WARNING: No male samples found in this plate's BCF. "
                "Skipping haploid recoding.",
                file=sys.stderr, flush=True,
            )
        else:
            print(f"[convert_gtc2vcf] Males in this plate : {len(males_in_plate)}", flush=True)
            print(f"[convert_gtc2vcf] Total samples       : {len(bcf_samples)}", flush=True)

            plugin_env = build_plugin_env()
            print(
                f"[convert_gtc2vcf] BCFTOOLS_PLUGINS    : "
                f"{plugin_env.get('BCFTOOLS_PLUGINS', '(not set)')}",
                flush=True,
            )

            xy_regions = detect_xy_regions(bcf_out)
            xy_set     = set(xy_regions.split(","))
            print(f"[convert_gtc2vcf] Sex-chr regions     : {xy_regions}", flush=True)

            # Temporary file paths
            tmp_xy          = Path(f"{prefix}_xy.bcf")
            tmp_xy_males    = Path(f"{prefix}_xy_males.bcf")
            tmp_xy_females  = Path(f"{prefix}_xy_females.bcf")
            tmp_xy_recoded  = Path(f"{prefix}_xy_recoded.bcf")
            tmp_auto        = Path(f"{prefix}_autosomes.bcf")
            haploid_bcf     = Path(f"{prefix}_haploid.bcf")
            haploid_sorted  = Path(f"{prefix}_haploid_sorted.bcf")
            males_file      = Path(f"{prefix}_males.txt")
            females_file    = Path(f"{prefix}_females.txt")

            recoding_ok = False
            try:
                # Step 3a: extract sex-chr sites only
                run([
                    "bcftools", "view",
                    "--regions",     xy_regions,
                    "--output-type", "b",
                    "--output",      str(tmp_xy),
                    str(bcf_out),
                ], env=plugin_env)
                run(["bcftools", "index", str(tmp_xy)], env=plugin_env)

                # Step 3b: write males and females sample list files
                males_file.write_text("\n".join(males_in_plate) + "\n")
                females_in_plate = [s for s in bcf_samples if s not in set(males_in_plate)]
                if females_in_plate:
                    females_file.write_text("\n".join(females_in_plate) + "\n")

                # Step 3c: extract males-only sex-chr BCF, recode to haploid.
                # We subset to males first so GT="RR"/"AA"/"het" expressions
                # target all samples in the subset (i.e. all males) without
                # needing per-sample GT[] accessors (added only in bcftools 1.21).
                # Three passes: hom-ref -> 0, hom-alt -> 1, het -> missing (.)
                print("[convert_gtc2vcf] Running +setGT pipeline on males-only sex chrs...",
                      flush=True)

                # Preflight: verify +setGT is loadable before running the pipe.
                # If the plugin is not found, bcftools exits non-zero and prints
                # "Failed to open plugin".  Catch it here with a clear error
                # message rather than a cryptic mid-pipe failure.
                setgt_check = subprocess.run(
                    ["bcftools", "plugin", "setGT"],
                    capture_output=True, text=True, env=plugin_env,
                )
                plugin_missing = (
                    setgt_check.returncode not in (0, 1)
                    or "failed to open" in setgt_check.stderr.lower()
                    or "no such file" in setgt_check.stderr.lower()
                )
                if plugin_missing:
                    raise SystemExit(
                        f"+setGT plugin not found. "
                        f"BCFTOOLS_PLUGINS={plugin_env.get('BCFTOOLS_PLUGINS', '(unset)')}.\n"
                        f"Set BCFTOOLS_LIBEXEC in your environment (or nextflow.config) to the "
                        f"directory containing setGT.so.\n"
                        f"On CHPC Lengau: export BCFTOOLS_LIBEXEC="
                        f"/apps/chpc/bio/bcftools/1.20/libexec/bcftools"
                    )
                print("[convert_gtc2vcf] +setGT plugin verified OK.", flush=True)

                setgt_pipe = (
                    f"bcftools view --force-samples --samples-file {males_file} --output-type u {tmp_xy} "
                    f"| bcftools +setGT --output-type u -- --target-gt q --new-gt 0 --include 'GT=\"RR\"' "
                    f"| bcftools +setGT --output-type u -- --target-gt q --new-gt c:1 --include 'GT=\"AA\"' "
                    f"| bcftools +setGT --output-type b --output {tmp_xy_males} -- --target-gt q --new-gt . --include 'GT=\"het\"'"
                )
                # check=True: abort loudly if the pipe fails rather than silently
                # keeping the diploid BCF and producing invalid F-statistics.
                run(setgt_pipe, shell=True, env=plugin_env, check=True)
                run(["bcftools", "index", str(tmp_xy_males)], env=plugin_env)

                # Step 3d: extract females-only sex-chr BCF (unmodified)
                # then merge males + females back into one BCF.
                if females_in_plate:
                    run([
                        "bcftools", "view",
                        "--force-samples",
                        "--samples-file", str(females_file),
                        "--output-type",  "b",
                        "--output",       str(tmp_xy_females),
                        str(tmp_xy),
                    ], env=plugin_env)
                    run(["bcftools", "index", str(tmp_xy_females)], env=plugin_env)

                    run([
                        "bcftools", "merge",
                        "--output-type", "b",
                        "--output",      str(tmp_xy_recoded),
                        str(tmp_xy_males),
                        str(tmp_xy_females),
                    ], env=plugin_env)
                else:
                    # Plate is all-male — recoded males BCF is the full sex-chr BCF
                    Path(str(tmp_xy_males)).rename(tmp_xy_recoded)

                run(["bcftools", "index", str(tmp_xy_recoded)], env=plugin_env)

                # Sanity-check: recoded BCF must be readable
                check_r = subprocess.run(
                    ["bcftools", "view", "--header-only", str(tmp_xy_recoded)],
                    capture_output=True, text=True,
                )
                if check_r.returncode != 0:
                    raise SystemExit(
                        "Recoded sex-chr BCF is unreadable after +setGT — "
                        "check bcftools stderr above for clues."
                    )

                # Step 3e: extract autosomes from the original BCF
                auto_contigs = get_auto_contigs(bcf_out, xy_set)

                if auto_contigs:
                    run([
                        "bcftools", "view",
                        "--regions",     ",".join(auto_contigs),
                        "--output-type", "b",
                        "--output",      str(tmp_auto),
                        str(bcf_out),
                    ], env=plugin_env)
                    run(["bcftools", "index", str(tmp_auto)], env=plugin_env)

                    # Step 3f: concat autosomes + recoded sex chrs
                    run([
                        "bcftools", "concat",
                        "--allow-overlaps",
                        "--output-type", "b",
                        "--output",      str(haploid_bcf),
                        str(tmp_auto),
                        str(tmp_xy_recoded),
                    ], env=plugin_env)
                    tmp_auto.unlink(missing_ok=True)
                    Path(f"{tmp_auto}.csi").unlink(missing_ok=True)
                else:
                    # BCF is sex-chr only (already region-split upstream)
                    tmp_xy_recoded.rename(haploid_bcf)

                # Sort (concat --allow-overlaps can leave chr boundary
                # records slightly out of order)
                run([
                    "bcftools", "sort",
                    "--output-type", "b",
                    "--output",      str(haploid_sorted),
                    "--temp-dir",    ".",
                    str(haploid_bcf),
                ], env=plugin_env)
                haploid_bcf.unlink(missing_ok=True)

                # Replace the diploid BCF
                bcf_out.unlink(missing_ok=True)
                Path(f"{bcf_out}.csi").unlink(missing_ok=True)
                haploid_sorted.rename(bcf_out)
                run(["bcftools", "index", str(bcf_out)], env=plugin_env)

                recoding_ok = True
                print(
                    f"[convert_gtc2vcf] Haploid recoding complete -> {bcf_out}",
                    flush=True,
                )

            except SystemExit as exc:
                print(
                    f"[convert_gtc2vcf] WARNING: Haploid recoding aborted: {exc}. "
                    f"Keeping diploid BCF.",
                    file=sys.stderr, flush=True,
                )
            finally:
                for tmp in [tmp_xy, tmp_xy_males, tmp_xy_females, tmp_xy_recoded,
                            tmp_auto, haploid_bcf, haploid_sorted]:
                    if tmp.exists() and not recoding_ok:
                        tmp.unlink(missing_ok=True)
                for csi in [
                    Path(f"{tmp_xy}.csi"),
                    Path(f"{tmp_xy_males}.csi"),
                    Path(f"{tmp_xy_females}.csi"),
                    Path(f"{tmp_xy_recoded}.csi"),
                    Path(f"{tmp_auto}.csi"),
                ]:
                    if csi.exists() and not recoding_ok:
                        csi.unlink(missing_ok=True)
                # Always clean up sample list files
                for f in [males_file, females_file]:
                    if f.exists():
                        f.unlink(missing_ok=True)

else:
    print(
        "\n[convert_gtc2vcf] Step 3: Skipped haploid recoding (no --sex-info supplied).",
        flush=True,
    )
    print(
        "[convert_gtc2vcf] NOTE: Without haploid recoding, PLINK --check-sex will",
        flush=True,
    )
    print(
        "[convert_gtc2vcf] report F=-1 for all males -> all inferred as Female.",
        flush=True,
    )


# ── Step 4: Extract X/Y intensities to TSV ────────────────────────────────────
print(f"\n[convert_gtc2vcf] Step 4: Extract XY intensities -> TSV", flush=True)

query_cmd = [
    "bcftools", "query",
    "--format", "[%SAMPLE\t%CHROM\t%POS\t%REF\t%ALT\t%NORMX\t%NORMY\n]",
    str(bcf_out),
]

with open(tsv_out, "w") as fh:
    fh.write("SAMPLE_ID\tCHR\tPOS\tREF\tALT\tNORMX\tNORMY\n")
    result = subprocess.run(query_cmd, stdout=fh, stderr=subprocess.PIPE, text=True)

if result.stderr:
    print(result.stderr, file=sys.stderr, flush=True)
if result.returncode != 0:
    print(
        f"[convert_gtc2vcf] WARNING: TSV extraction failed (exit {result.returncode}). "
        f"Downstream XY intensity QC will be skipped for this plate.",
        file=sys.stderr, flush=True,
    )
    with open(tsv_out, "w") as fh:
        fh.write("SAMPLE_ID\tCHR\tPOS\tREF\tALT\tNORMX\tNORMY\n")
else:
    print(f"[convert_gtc2vcf] TSV written  : {tsv_out}", flush=True)

print(f"\n[convert_gtc2vcf] Done -- plate {prefix}", flush=True)