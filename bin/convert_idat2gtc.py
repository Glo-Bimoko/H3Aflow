"""
convert_idat2gtc.py
===================
Converts a single sample's idat files to a GTC file using bcftools +idat2gtc.

bcftools +idat2gtc --idats <dir> writes output GTC files into an output
directory named automatically as <barcode>_<position>.gtc. This script
runs the conversion then renames the output to the desired sample name.

If the output GTC file already exists at the destination AND passes a health
check, the conversion is skipped and the existing file is used as-is.
If the file exists but is broken, it is deleted and reconverted from the
original idats.

GTC health check
----------------
check_gtc_health() runs all four checks on every call and returns a structured
GtcHealthReport.  All checks are always executed so the log always shows the
complete picture even when an earlier check has already failed.

Three checks are *mandatory* — a failure in any one marks the file as broken
and triggers deletion + reconversion:
  1. File size       — must be >= GTC_MIN_BYTES (1 MB).  A valid GTC for a
                       standard Illumina array is several MB; anything smaller
                       is a truncated write or a zero-byte stub left by a
                       crashed job.
  2. Magic bytes     — the first 3 bytes must be the ASCII string "gtc"
                       (0x67 0x74 0x63), the Illumina GTC format identifier.
  3. GTC header parse — the file is opened and its binary Table-of-Contents
                       (TOC) header is parsed natively in Python, following
                       the official Illumina GTC format specification:
                         bytes 0-2  : "gtc" magic
                         byte  3    : version (int8)
                         bytes 4-7  : number of TOC entries (int32 LE)
                         N × 6 bytes: TOC entry = int16 id + uint32 value/offset
                       TOC entry id=1 holds num_snps directly as its uint32
                       value field (not a file offset).  The check verifies:
                         a) the header is readable and well-formed
                         b) num_snps > 0 (file is not empty/zeroed-out)
                         c) the actual file size is >= the minimum size implied
                            by num_snps, computed as:
                              header_size + num_snps × bytes_per_snp_floor
                            where bytes_per_snp_floor covers genotypes (1 byte)
                            + raw X (2 bytes) + raw Y (2 bytes) = 5 bytes/SNP
                            minimum, plus a modest fixed header allowance.
                       This is a pure Python, zero-subprocess check that reads
                       only the first ~few hundred bytes of the file and is
                       therefore very fast even for large GTC files.

One check is *informational* — its result is always logged but never triggers
reconversion on its own:
  4. Version byte    — byte index 3 should be 3-5 (Illumina GTC versions in
                       active use).  Values outside this range may indicate an
                       unusual tool version but do not necessarily mean the
                       file is unreadable.

GTC binary layout reference:
    Illumina BeadArrayFiles — GenotypeCalls.py (official Illumina parser)
    https://github.com/Illumina/BeadArrayFiles

Usage:
    python convert_idat2gtc.py \\
        --bpm       /path/to/manifest.bpm \\
        --egt       /path/to/clusters.egt \\
        --idats     /path/to/per_sample_idat_dir \\
        --output    sample_id.gtc \\
        [--gtc-dir  /path/to/published/gtc/dir]
"""

import argparse
import shutil
import struct
import subprocess
import sys
from dataclasses import dataclass, field
from pathlib import Path

# ── Constants ──────────────────────────────────────────────────────────────────
GTC_MAGIC       = b"gtc"           # first 3 bytes of every valid GTC file
GTC_MIN_VERSION = 3                # earliest known GTC version byte
GTC_MAX_VERSION = 5                # latest known GTC version byte
GTC_MIN_BYTES   = 1 * 1024 * 1024  # 1 MB — below this the file is considered truncated

# TOC entry id for num_snps (value stored directly in the offset field of the TOC entry)
_GTC_TOC_ID_NUM_SNPS = 1

# Minimum bytes per SNP used to estimate the smallest possible valid file size.
# Covers: genotype byte (1) + raw X uint16 (2) + raw Y uint16 (2) = 5 bytes,
# plus the 4-byte count prefix for each array.  This is a conservative floor;
# real GTCs are much larger (they also carry scores, BAF, LRR, transforms, etc).
_GTC_BYTES_PER_SNP_FLOOR = 5
# Fixed overhead for the GTC header itself (magic + version + toc count + toc
# entries + string fields).  500 bytes is a very conservative lower bound.
_GTC_HEADER_OVERHEAD = 500


# ── Health-check result dataclasses ───────────────────────────────────────────
@dataclass
class CheckResult:
    name: str
    passed: bool
    mandatory: bool
    detail: str


@dataclass
class GtcHealthReport:
    path: Path
    checks: list = field(default_factory=list)  # list[CheckResult]

    @property
    def is_healthy(self):
        """True only when every mandatory check passed."""
        return all(c.passed for c in self.checks if c.mandatory)

    @property
    def failed_mandatory(self):
        return [c for c in self.checks if c.mandatory and not c.passed]

    @property
    def failed_informational(self):
        return [c for c in self.checks if not c.mandatory and not c.passed]

    def log(self, label=""):
        """Print a structured summary of all four checks."""
        prefix = "[convert_idat2gtc][health-check]"
        tag = (prefix + " " + label) if label else prefix
        print(tag + " File: " + str(self.path), flush=True)
        for c in self.checks:
            kind  = "MANDATORY    " if c.mandatory else "informational"
            state = "PASS" if c.passed else ("FAIL" if c.mandatory else "WARN")
            print("  [" + state + "] (" + kind + ") " + c.name + ": " + c.detail, flush=True)
        overall = "HEALTHY" if self.is_healthy else "BROKEN"
        print("  -> Overall: " + overall, flush=True)


# ── GTC header parser (pure Python, no subprocesses) ──────────────────────────
def _parse_gtc_header(path):
    """
    Parse the GTC binary TOC header and return (num_snps, version, toc_table)
    or raise an exception describing what went wrong.

    GTC binary layout (from official Illumina BeadArrayFiles spec):
      bytes 0-2  : "gtc" magic string
      byte  3    : file version (uint8)
      bytes 4-7  : number of TOC entries (int32 LE)
      N × 6 bytes: TOC entries, each = struct "<hI"
                     int16  id     : field identifier
                     uint32 value  : for ids 1/2/3 this IS the field value;
                                     for all other ids this is a file offset
    """
    with path.open("rb") as fh:
        # magic + version already validated by check 2, but re-read cleanly here
        header_start = fh.read(8)
        if len(header_start) < 8:
            raise ValueError("file too short to contain a GTC header (" + str(len(header_start)) + " bytes)")

        version = header_start[3]
        num_toc = struct.unpack_from("<I", header_start, 4)[0]

        if num_toc == 0 or num_toc > 10000:
            raise ValueError("implausible TOC entry count: " + str(num_toc) +
                             " (expected 1-10000); file may be zeroed or corrupt")

        toc_bytes = fh.read(num_toc * 6)
        if len(toc_bytes) < num_toc * 6:
            raise ValueError("TOC truncated: expected " + str(num_toc * 6) +
                             " bytes for " + str(num_toc) + " entries, got " + str(len(toc_bytes)))

        toc_table = {}
        for i in range(num_toc):
            entry_id, value = struct.unpack_from("<hI", toc_bytes, i * 6)
            toc_table[entry_id] = value

        if _GTC_TOC_ID_NUM_SNPS not in toc_table:
            raise ValueError("TOC does not contain a num_snps entry (id=" +
                             str(_GTC_TOC_ID_NUM_SNPS) + "); file may be corrupt")

        num_snps = toc_table[_GTC_TOC_ID_NUM_SNPS]
        return num_snps, version, toc_table


# ── Four-check health function ─────────────────────────────────────────────────
def check_gtc_health(path):
    """
    Run all four GTC health checks against *path* and return a GtcHealthReport.

    All checks are always executed so the log always shows the complete picture,
    even when an earlier check has already failed.

    Parameters
    ----------
    path : Path   GTC file to inspect.

    Note: the bcftools read test was replaced by a native Python GTC header
    parse (check 3) which reads only the first ~few hundred bytes of the file
    and is therefore much cheaper than a full bcftools conversion.
    """
    report = GtcHealthReport(path=path)

    # ── Pre-flight: file must exist and be a regular file ─────────────────────
    if not path.exists() or not path.is_file():
        skipped = "skipped — file does not exist or is not a regular file"
        for name, mandatory in [
            ("File size",        True),
            ("Magic bytes",      True),
            ("GTC header parse", True),
            ("Version byte",     False),
        ]:
            report.checks.append(CheckResult(
                name=name, passed=False, mandatory=mandatory, detail=skipped,
            ))
        return report

    # ── Check 1 (mandatory): File size ────────────────────────────────────────
    size = path.stat().st_size
    if size >= GTC_MIN_BYTES:
        c1 = CheckResult(
            name="File size", passed=True, mandatory=True,
            detail=str(size) + " bytes (>= " + str(GTC_MIN_BYTES) + " byte minimum)",
        )
    else:
        c1 = CheckResult(
            name="File size", passed=False, mandatory=True,
            detail=(str(size) + " bytes < " + str(GTC_MIN_BYTES) + " byte minimum; "
                    "likely a truncated write or zero-byte placeholder"),
        )
    report.checks.append(c1)

    # ── Read header bytes (shared by checks 2 & 4) ───────────────────────────
    header = b""
    header_error = None
    try:
        with path.open("rb") as fh:
            header = fh.read(4)
    except OSError as exc:
        header_error = str(exc)

    # ── Check 2 (mandatory): Magic bytes ──────────────────────────────────────
    if header_error:
        c2 = CheckResult(
            name="Magic bytes", passed=False, mandatory=True,
            detail="cannot read file header: " + header_error,
        )
    elif len(header) < 3:
        c2 = CheckResult(
            name="Magic bytes", passed=False, mandatory=True,
            detail="header too short (" + str(len(header)) + " bytes); file severely truncated",
        )
    elif header[:3] == GTC_MAGIC:
        c2 = CheckResult(
            name="Magic bytes", passed=True, mandatory=True,
            detail="found expected magic b'gtc' (0x" + GTC_MAGIC.hex() + ")",
        )
    else:
        actual = header[:3]
        c2 = CheckResult(
            name="Magic bytes", passed=False, mandatory=True,
            detail=("expected b'gtc' (0x" + GTC_MAGIC.hex() + "), "
                    "got " + repr(actual) + " (0x" + actual.hex() + "); "
                    "file is not a GTC or has been overwritten"),
        )
    report.checks.append(c2)

    # ── Check 3 (mandatory): GTC header parse ─────────────────────────────────
    # Parse the binary TOC header natively in Python — reads only the first
    # ~few hundred bytes.  Verifies:
    #   (a) the header is well-formed and parseable
    #   (b) num_snps > 0
    #   (c) actual file size >= minimum size implied by num_snps
    # This replaces the previous bcftools read test, which paid full conversion
    # cost just to validate an existing file.
    try:
        num_snps, _version, _toc = _parse_gtc_header(path)

        if num_snps == 0:
            c3 = CheckResult(
                name="GTC header parse", passed=False, mandatory=True,
                detail="num_snps is 0 in TOC header; file is empty or was zeroed out",
            )
        else:
            min_expected = _GTC_HEADER_OVERHEAD + num_snps * _GTC_BYTES_PER_SNP_FLOOR
            if size < min_expected:
                c3 = CheckResult(
                    name="GTC header parse", passed=False, mandatory=True,
                    detail=("header claims num_snps=" + str(num_snps) +
                            " but file size " + str(size) + " bytes is below the minimum "
                            "expected size of " + str(min_expected) + " bytes "
                            "(" + str(num_snps) + " SNPs × " + str(_GTC_BYTES_PER_SNP_FLOOR) +
                            " bytes/SNP + " + str(_GTC_HEADER_OVERHEAD) + " byte header overhead); "
                            "file was truncated mid-write"),
                )
            else:
                c3 = CheckResult(
                    name="GTC header parse", passed=True, mandatory=True,
                    detail=("TOC header parsed successfully; num_snps=" + str(num_snps) +
                            ", file size " + str(size) + " bytes is consistent "
                            "(minimum expected " + str(min_expected) + " bytes)"),
                )
    except Exception as exc:
        c3 = CheckResult(
            name="GTC header parse", passed=False, mandatory=True,
            detail="failed to parse GTC TOC header: " + str(exc),
        )
    report.checks.append(c3)

    # ── Check 4 (informational): Version byte ─────────────────────────────────
    if header_error or len(header) < 4:
        c4 = CheckResult(
            name="Version byte", passed=False, mandatory=False,
            detail="skipped — header unreadable (see check 2)",
        )
    else:
        version = header[3]
        if GTC_MIN_VERSION <= version <= GTC_MAX_VERSION:
            c4 = CheckResult(
                name="Version byte", passed=True, mandatory=False,
                detail=("version " + str(version) +
                        " (supported range " + str(GTC_MIN_VERSION) + "-" + str(GTC_MAX_VERSION) + ")"),
            )
        else:
            c4 = CheckResult(
                name="Version byte", passed=False, mandatory=False,
                detail=("version byte is " + str(version) + ", outside expected range "
                        + str(GTC_MIN_VERSION) + "-" + str(GTC_MAX_VERSION) + "; "
                        "may be an unusual tool version — treated as informational only"),
            )
    report.checks.append(c4)

    return report


# ── Bulk scan function ────────────────────────────────────────────────────────
def bulk_scan_gtc_dir(gtc_dir):
    """
    Scan every *.gtc file in gtc_dir, delete any that fail the health check,
    and print a summary to stdout.

    This runs exactly once per pipeline invocation — whichever IDAT_TO_GTC task
    wins the atomic .scan_done marker creation runs this; all others skip it.

    Parameters
    ----------
    gtc_dir : Path   The published GTC output directory (e.g. results/gtc/).
    """
    gtc_files = sorted(gtc_dir.glob("*.gtc"))

    if not gtc_files:
        print(
            "[convert_idat2gtc][bulk-scan] No GTC files found in " + str(gtc_dir) + " — nothing to scan.",
            flush=True,
        )
        return

    print(
        "[convert_idat2gtc][bulk-scan] Scanning " + str(len(gtc_files)) +
        " GTC file(s) in " + str(gtc_dir) + " ...",
        flush=True,
    )

    healthy_count = 0
    broken_count  = 0
    deleted_count = 0

    for gtc_path in gtc_files:
        report = check_gtc_health(gtc_path)
        if report.is_healthy:
            healthy_count += 1
        else:
            broken_count += 1
            failed_names = ", ".join(c.name for c in report.failed_mandatory)
            print(
                "[convert_idat2gtc][bulk-scan] BROKEN: " + gtc_path.name +
                " (failed: " + failed_names + ") — deleting.",
                flush=True,
            )
            # Log each failing check for full visibility
            for c in report.failed_mandatory:
                print(
                    "  [FAIL] " + c.name + ": " + c.detail,
                    flush=True,
                )
            try:
                gtc_path.unlink()
                deleted_count += 1
                print(
                    "[convert_idat2gtc][bulk-scan] Deleted: " + str(gtc_path),
                    flush=True,
                )
            except OSError as exc:
                print(
                    "[convert_idat2gtc][bulk-scan] WARNING: Could not delete " +
                    str(gtc_path) + ": " + str(exc),
                    flush=True,
                )

    print(
        "[convert_idat2gtc][bulk-scan] Complete — "
        + str(healthy_count) + " healthy, "
        + str(broken_count)  + " broken ("
        + str(deleted_count) + " deleted).",
        flush=True,
    )


# ── CLI ────────────────────────────────────────────────────────────────────────
parser = argparse.ArgumentParser()
parser.add_argument("--bpm",       required=True)
parser.add_argument("--egt",       required=True)
parser.add_argument("--idats",     required=True)
parser.add_argument("--output",    required=True)
parser.add_argument(
    "--gtc-dir",
    default=None,
    help=(
        "Published GTC output directory (e.g. params.outdir/gtc). "
        "When supplied, the script checks this directory for an existing "
        "<sample_id>.gtc before running conversion."
    ),
)
args = parser.parse_args()

bpm    = Path(args.bpm).resolve()
egt    = Path(args.egt).resolve()
idats  = Path(args.idats).resolve()
output = Path(args.output).resolve()   # desired final path e.g. SAMPLE_ID.gtc

for f, label in [(bpm, "BPM"), (egt, "EGT"), (idats, "idat directory")]:
    if not f.exists():
        sys.exit("[convert_idat2gtc] ERROR: " + label + " not found: " + str(f))

# ── Bulk scan: runs once across all parallel tasks ────────────────────────────
# The first task to atomically create .scan_done wins and runs the bulk scan.
# All other tasks see the file already exists and skip straight to their own
# per-sample check below.
# open(..., 'x') is an atomic exclusive-create on POSIX filesystems — only one
# process succeeds even if many attempt it simultaneously.
if args.gtc_dir:
    gtc_dir_path   = Path(args.gtc_dir).resolve()
    scan_done_flag = gtc_dir_path / ".scan_done"
    gtc_dir_path.mkdir(parents=True, exist_ok=True)
    try:
        with open(scan_done_flag, 'x') as _flag_fh:
            pass  # we won the race — run the bulk scan
        print(
            "[convert_idat2gtc][bulk-scan] This task won the scan lock — "
            "running bulk scan before per-sample conversion.",
            flush=True,
        )
        bulk_scan_gtc_dir(gtc_dir_path)
    except FileExistsError:
        print(
            "[convert_idat2gtc][bulk-scan] Scan already done (or in progress) — skipping.",
            flush=True,
        )

# ── Skip-or-reconvert logic ───────────────────────────────────────────────────
# Priority order:
#   1. Published GTC directory (permanent results location, checked first).
#   2. Local output path (Nextflow work directory for this task).
#
# For whichever candidate is found first, run the full four-check health report.
#   All mandatory checks pass → skip conversion; copy to output if needed.
#   Any mandatory check fails → delete broken file; fall through to conversion.
#   Not found                 → fall through to conversion.
sample_id = output.stem

existing_gtc = None

if args.gtc_dir:
    candidate = Path(args.gtc_dir).resolve() / (sample_id + ".gtc")
    if candidate.exists():
        existing_gtc = candidate

if existing_gtc is None and output.exists():
    existing_gtc = output

if existing_gtc is not None:
    report = check_gtc_health(existing_gtc)
    report.log(label=sample_id)

    # Always surface informational warnings even for an otherwise healthy file.
    for c in report.failed_informational:
        print(
            "[convert_idat2gtc] WARNING (" + sample_id + "): "
            + c.name + " — " + c.detail,
            file=sys.stderr, flush=True,
        )

    if report.is_healthy:
        print(
            "[convert_idat2gtc] SKIP  " + sample_id + ": GTC exists and is healthy "
            "(" + str(existing_gtc) + ")",
            flush=True,
        )
        # Copy to the expected output path if the healthy file lives elsewhere
        # (e.g. in the published results dir rather than the Nextflow work dir).
        if existing_gtc != output:
            shutil.copy2(str(existing_gtc), str(output))
            print(
                "[convert_idat2gtc] Copied existing GTC -> " + str(output),
                flush=True,
            )
        sys.exit(0)

    else:
        # One or more mandatory checks failed — delete the broken file and reconvert.
        failed_names = ", ".join(c.name for c in report.failed_mandatory)
        print(
            "[convert_idat2gtc] WARNING: Existing GTC for " + sample_id + " is broken "
            "(failed mandatory checks: " + failed_names + ").\n"
            "  Deleting broken file: " + str(existing_gtc),
            file=sys.stderr, flush=True,
        )
        try:
            existing_gtc.unlink()
            print(
                "[convert_idat2gtc] Deleted broken GTC: " + str(existing_gtc),
                flush=True,
            )
        except OSError as exc:
            # Non-fatal: warn and continue; conversion will overwrite.
            print(
                "[convert_idat2gtc] WARNING: Could not delete broken GTC "
                "(" + str(existing_gtc) + "): " + str(exc) + ". "
                "Proceeding with conversion anyway.",
                file=sys.stderr, flush=True,
            )
        # Also remove the local output copy if it is a different path and exists.
        if existing_gtc != output and output.exists():
            try:
                output.unlink()
            except OSError:
                pass

        # Fall through to conversion below.

# ── Conversion (with retry) ───────────────────────────────────────────────────
# bcftools +idat2gtc writes GTC files into a directory.
# Up to MAX_ATTEMPTS attempts are made.  If a newly produced GTC fails the
# health check the temp directory is cleaned up and the conversion is retried.
# After all attempts are exhausted the sample is abandoned with a WARNING so
# the rest of the pipeline can continue.
MAX_ATTEMPTS = 3

for attempt in range(1, MAX_ATTEMPTS + 1):
    attempt_tag = "[attempt " + str(attempt) + "/" + str(MAX_ATTEMPTS) + "]"

    gtc_outdir = output.parent() / ("_gtc_tmp_" + output.stem)  # always absolute — anchored to Nextflow work dir
    shutil.rmtree(str(gtc_outdir), ignore_errors=True)  # clean any leftover from prior attempt
    gtc_outdir.mkdir(parents=True, exist_ok=True)

    cmd = [
        "bcftools", "+idat2gtc",
        "--bpm",    str(bpm),
        "--egt",    str(egt),
        "--idats",  str(idats),
        "--output", str(gtc_outdir),
    ]

    print("[convert_idat2gtc] " + attempt_tag + " idat dir : " + str(idats), flush=True)
    print("[convert_idat2gtc] " + attempt_tag + " gtc outdir: " + str(gtc_outdir), flush=True)
    print("[convert_idat2gtc] " + attempt_tag + " Running  : " + " ".join(cmd) + "\n", flush=True)

    result = subprocess.run(cmd, capture_output=True, text=True)

    if result.stdout:
        print(result.stdout, flush=True)
    if result.stderr:
        print(result.stderr, file=sys.stderr, flush=True)

    if result.returncode != 0:
        # Check if the failure was due to a missing idat file or numerical failure.
        stderr_text = result.stderr if result.stderr else ""
        if ("No such file or directory" in stderr_text
                or "Could not open" in stderr_text
                or "Error while running linsolve" in stderr_text
                or "linsolve" in stderr_text.lower()):
            print(
                "[convert_idat2gtc] WARNING: Some samples could not be processed in "
                + str(idats) + ". "
                "Reason may be missing idat files or numerical failure during normalization. "
                "Skipping affected samples and continuing with available pairs.",
                file=sys.stderr, flush=True
            )
            # Check if any GTCs were produced before the failure.
            gtc_files = list(gtc_outdir.glob("*.gtc")) if gtc_outdir.exists() else []
            if not gtc_files:
                shutil.rmtree(str(gtc_outdir), ignore_errors=True)
                sys.exit(
                    "[convert_idat2gtc] ERROR: No GTC files produced and idat files missing. "
                    "This sample group cannot be processed."
                )
            print(
                "[convert_idat2gtc] Partial success: " + str(len(gtc_files)) + " GTC(s) produced.",
                flush=True
            )
        else:
            sys.exit(
                "[convert_idat2gtc] ERROR: bcftools +idat2gtc failed "
                "(exit code " + str(result.returncode) + ")"
            )

    # Find the GTC file written into the temp output directory.
    gtc_files = list(gtc_outdir.glob("*.gtc"))
    if not gtc_files:
        sys.exit(
            "[convert_idat2gtc] ERROR: No GTC file produced in " + str(gtc_outdir) + "\n"
            "  Contents: " + str(list(gtc_outdir.iterdir()))
        )

    if len(gtc_files) > 1:
        print(
            "[convert_idat2gtc] WARNING: Multiple GTC files found: " + str(gtc_files) + "\n"
            "  Using first: " + str(gtc_files[0]),
            file=sys.stderr
        )

    produced_gtc = gtc_files[0]
    print("[convert_idat2gtc] " + attempt_tag + " Produced  : " + produced_gtc.name, flush=True)

    # ── Validate the freshly produced GTC ─────────────────────────────────────
    fresh_report = check_gtc_health(produced_gtc)
    fresh_report.log(label=sample_id + " " + attempt_tag)

    for c in fresh_report.failed_informational:
        print(
            "[convert_idat2gtc] WARNING (" + sample_id + "): "
            + c.name + " — " + c.detail,
            file=sys.stderr, flush=True,
        )

    if fresh_report.is_healthy:
        # Good GTC — move to the final output path and exit the retry loop.
        shutil.move(str(produced_gtc), str(output))
        shutil.rmtree(str(gtc_outdir), ignore_errors=True)
        break

    # Health check failed — clean up and decide whether to retry or abandon.
    failed_names = ", ".join(c.name for c in fresh_report.failed_mandatory)
    shutil.rmtree(str(gtc_outdir), ignore_errors=True)

    if attempt < MAX_ATTEMPTS:
        print(
            "[convert_idat2gtc] WARNING " + attempt_tag + ": Produced GTC for " + sample_id
            + " failed health check (failed mandatory checks: " + failed_names + "). "
            "Retrying conversion.",
            file=sys.stderr, flush=True,
        )
    else:
        # All attempts exhausted — abandon this sample with a warning.
        print(
            "[convert_idat2gtc] WARNING: " + sample_id + " ABANDONED after " + str(MAX_ATTEMPTS)
            + " conversion attempt(s). The produced GTC failed health checks every time "
            "(failed mandatory checks: " + failed_names + "). "
            "This sample will be missing from downstream results. "
            "Source idats: " + str(idats),
            file=sys.stderr, flush=True,
        )
        sys.exit(0)  # exit 0 so Nextflow does not retry the task as a failure

else:
    # Loop completed without break — should be unreachable, but guard anyway.
    sys.exit("[convert_idat2gtc] ERROR: Retry loop exited unexpectedly for " + sample_id)

if not output.exists():
    sys.exit("[convert_idat2gtc] ERROR: Output GTC not found after rename: " + str(output))

print("[convert_idat2gtc] Done -> " + str(output), flush=True)