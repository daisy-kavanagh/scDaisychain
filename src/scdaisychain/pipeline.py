#!/usr/bin/env python3
"""Full XCI pipeline runner."""

from __future__ import annotations

import argparse
import shlex
import shutil
import subprocess
import sys
from pathlib import Path
from typing import Iterable, List


MODULES = {
    "count": "scdaisychain.count_snps",
    "phase": "scdaisychain.phase_x",
    "tag": "scdaisychain.tag_bam",
    "split": "scdaisychain.split_bam",
    "matrix": "scdaisychain.make_matrices",
}


def lc_tag(lower_cutoff: float) -> str:
    """Keep naming consistent with phase_x.py."""
    return f"lc{float(lower_cutoff):.2f}"


def shell_join(cmd: Iterable[object]) -> str:
    return " ".join(shlex.quote(str(x)) for x in cmd)


def ensure_file(path: Path, label: str) -> None:
    if not path.exists():
        raise FileNotFoundError(f"{label} does not exist: {path}")
    if not path.is_file():
        raise FileNotFoundError(f"{label} is not a file: {path}")


def ensure_dir(path: Path, label: str) -> None:
    if not path.exists():
        raise FileNotFoundError(f"{label} does not exist: {path}")
    if not path.is_dir():
        raise FileNotFoundError(f"{label} is not a directory: {path}")


def add_flag(cmd: List[str], flag: str, enabled: bool) -> None:
    if enabled:
        cmd.append(flag)


def add_opt(cmd: List[str], flag: str, value) -> None:
    if value is not None:
        cmd.extend([flag, str(value)])


def module_cmd(python_executable: str, module_key: str) -> List[str]:
    return [python_executable, "-m", MODULES[module_key]]


def run_step(name: str, cmd: List[str], commands_log: Path, dry_run: bool = False) -> None:
    print(f"\n=== {name} ===", flush=True)
    print(shell_join(cmd), flush=True)

    with commands_log.open("a") as handle:
        handle.write(f"\n# {name}\n")
        handle.write(shell_join(cmd) + "\n")

    if dry_run:
        return

    subprocess.run(cmd, check=True)


def maybe_skip(step_name: str, expected_outputs: Iterable[Path], force: bool) -> bool:
    expected = list(expected_outputs)
    if expected and all(p.exists() for p in expected) and not force:
        print(f"\n=== {step_name} ===")
        print("Skipping because expected output already exists:")
        for p in expected:
            print(f"  {p}")
        print("Use --force to rerun this step.")
        return True
    return False


def add_run_arguments(p: argparse.ArgumentParser) -> argparse.ArgumentParser:
    # Core inputs/outputs
    p.add_argument("--bam", required=True, help="Input coordinate-sorted/indexed BAM/CRAM used for counting and tagging.")
    p.add_argument("--vcf", required=True, help="Input VCF/BCF of variants for SNP counting. Must be indexed for region fetch.")
    p.add_argument("--gtf", required=True, help="GTF used for count annotation/filtering and chrX matrix features.")
    p.add_argument("--outdir", required=True, help="Top-level output directory for this run.")
    p.add_argument("--run-name", default=None, help="Output prefix. Default: input BAM stem.")
    p.add_argument("--sample-name", default=None, help="Sample name written to the phased VCF. Default: --run-name.")
    p.add_argument("--python", default=sys.executable, help="Python executable used to run helper modules. Default: current Python.")

    # Resume/skip controls
    p.add_argument("--skip-count", action="store_true")
    p.add_argument("--skip-phase", action="store_true")
    p.add_argument("--skip-tag", action="store_true")
    p.add_argument("--skip-split", action="store_true")
    p.add_argument("--skip-matrices", action="store_true")
    p.add_argument("--force", action="store_true", help="Rerun steps even if expected outputs exist.")
    p.add_argument("--dry-run", action="store_true", help="Print commands without executing them.")

    # Manually supplied intermediate paths for resuming
    p.add_argument("--counts-tsv", default=None, help="Existing or desired SNP count TSV path.")
    p.add_argument("--active-x-csv", default=None, help="Existing or desired active-X CSV path.")
    p.add_argument("--phased-vcf", default=None, help="Existing or desired phased VCF.gz path for tagging.")
    p.add_argument("--tagged-bam", default=None, help="Existing or desired tagged BAM path.")
    p.add_argument("--split-bam-dir", default=None, help="Existing or desired split BAM directory.")
    p.add_argument("--matrices-outdir", default=None, help="Existing or desired matrix output directory.")

    # SNP counting options
    group = p.add_argument_group("SNP counting options")
    group.add_argument("--chrom", default="chrX", help="Chromosome to count. Default: chrX. Ignored if --region is set.")
    group.add_argument("--region", default=None, help="Region to count, e.g. chrX:1-20000000.")
    group.add_argument("--count-sample", default=None, help="VCF sample used by --het-only in the counting script.")
    group.add_argument("--het-only", action="store_true", help="Count only heterozygous SNPs from the input VCF.")
    group.add_argument("--procs", type=int, default=1, help="Parallel region workers for counting. Default: 1.")
    group.add_argument("--threads", type=int, default=1, help="HTSlib threads per counting worker. Default: 1.")
    group.add_argument("--window", type=int, default=5_000_000, help="Window size for parallel counting. Default: 5 Mb.")
    group.add_argument("--cb-tag", default="CB", help="Cell barcode tag. Default: CB.")
    group.add_argument("--min-mapq", type=int, default=20)
    group.add_argument("--min-bq", type=int, default=10)
    group.add_argument("--max-depth", type=int, default=100000)
    group.add_argument("--allow-secondary", action="store_true")
    group.add_argument("--count-duplicates", action="store_true")
    group.add_argument("--skip-indels", action="store_true")
    group.add_argument("--ignore-overlaps", action="store_true")
    group.add_argument("--slop", type=int, default=1000)
    group.add_argument("--max-nearest-dist", type=int, default=10000)
    group.add_argument("--drop-conflicts", action="store_true")
    group.add_argument("--drop-multi-tsv-and-gtf", action="store_true")
    group.add_argument("--no-dropped-read-outputs", action="store_true", help="Do not write dropped-read evidence files from the count step.")
    group.add_argument("--read-evidence-out", default=None, help="Optional all-read SNP evidence TSV from the count step.")
    group.add_argument("--count-debug", action="store_true")

    # Phasing options
    group = p.add_argument_group("Phasing options")
    group.add_argument("--min-reads", type=int, default=10, help="Minimum total reads per SNP for phasing. Default: 20.")
    group.add_argument("--lower-cutoff", type=float, default=0.01, help="Lower ref-fraction cutoff. Default: 0.10.")
    group.add_argument("--upper-cutoff", type=float, default=None)
    group.add_argument("--partition-mode", choices=["unweighted", "weighted", "expression"], default="weighted")
    group.add_argument("--edge-weight-mode", choices=["binary", "cosine", "overlap"], default=None)
    group.add_argument("--discordant-tsv", default=None)
    group.add_argument("--stage2-mode", choices=["cell", "expression"], default="cell")
    group.add_argument("--stage2-min-evidence", type=float, default=1.0)
    group.add_argument("--stage2-tie-action", choices=["random", "skip"], default="random")

    # BAM tagging/splitting options
    group = p.add_argument_group("BAM tagging and splitting options")
    group.add_argument("--tag-mode", choices=["quality", "count"], default="quality")
    group.add_argument("--tag-min-base-qual", type=int, default=5)
    group.add_argument("--tag-min-score-diff", type=int, default=10)
    group.add_argument("--tag-min-total-score", type=int, default=10)
    group.add_argument("--tag-min-count-diff", type=int, default=1)
    group.add_argument("--tag-min-total-count", type=int, default=1)
    group.add_argument("--split-prefix", default=None, help="Prefix for split BAMs. Default: --run-name.")

    # Matrix options
    group = p.add_argument_group("Matrix options")
    group.add_argument("--original-gene-matrix-dir", default=None, help="Original gene matrix dir containing barcodes.tsv.gz/features.tsv.gz.")
    group.add_argument("--original-transcript-matrix-dir", default=None, help="Original transcript matrix dir containing barcodes.tsv.gz/features.tsv.gz.")
    group.add_argument("--gn-tag", default="GN")
    group.add_argument("--tr-tag", default="TR")
    group.add_argument("--transcript-id-type", choices=["auto", "transcript_id", "transcript_name"], default="auto")
    group.add_argument("--matrix-bam-pattern", default="*.bam")
    group.add_argument(
        "--matrix-exclude-dropped-reads",
        action=argparse.BooleanOptionalAction,
        default=True,
        help="Pass unique dropped read IDs from count filtering to matrix generation when available. Default: true.",
    )
    group.add_argument("--matrix-exclude-duplicates", action="store_true")
    group.add_argument("--matrix-keep-secondary", action="store_true")
    group.add_argument("--matrix-keep-supplementary", action="store_true")
    group.add_argument("--skip-transcript", action="store_true")

    return p


def parse_run_args(argv=None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Run the full scDaisychain/XCI pipeline.")
    add_run_arguments(parser)
    return parser.parse_args(argv)


def run_pipeline(args: argparse.Namespace) -> int:
    bam = Path(args.bam).expanduser().resolve()
    vcf = Path(args.vcf).expanduser().resolve()
    gtf = Path(args.gtf).expanduser().resolve()
    outdir = Path(args.outdir).expanduser().resolve()

    ensure_file(bam, "Input BAM/CRAM")
    ensure_file(vcf, "Input VCF/BCF")
    ensure_file(gtf, "Input GTF")

    run_name = args.run_name or bam.name.removesuffix(".bam").removesuffix(".cram")
    sample_name = args.sample_name or run_name
    split_prefix = args.split_prefix or run_name
    lctag = lc_tag(args.lower_cutoff)

    counts_dir = outdir / "01_counts"
    phase_dir = outdir / "02_phase"
    tagged_dir = outdir / "03_tagged_bam"
    split_dir = Path(args.split_bam_dir).expanduser().resolve() if args.split_bam_dir else outdir / "04_split_bams"
    matrix_dir = Path(args.matrices_outdir).expanduser().resolve() if args.matrices_outdir else outdir / "05_matrices"
    logs_dir = outdir / "logs"

    for d in [counts_dir, phase_dir, tagged_dir, split_dir, matrix_dir, logs_dir]:
        d.mkdir(parents=True, exist_ok=True)

    commands_log = logs_dir / "commands.sh"
    header = "#!/usr/bin/env bash\nset -euo pipefail\n"
    if args.dry_run:
        header = "#!/usr/bin/env bash\n# dry-run commands\nset -euo pipefail\n"
    commands_log.write_text(header)

    counts_tsv = Path(args.counts_tsv).expanduser().resolve() if args.counts_tsv else counts_dir / f"{run_name}.snp_counts.tsv"
    count_summary = counts_dir / f"{run_name}.snp_counts.summary.tsv"
    dropped_read_evidence = counts_dir / f"{run_name}.dropped_read_evidence.tsv"
    unique_dropped_read_ids = counts_dir / f"{run_name}.unique_dropped_read_ids.txt"
    dropped_read_summary = counts_dir / f"{run_name}.dropped_read_summary.tsv"

    default_active_x_csv = phase_dir / f"haplotype_sums_df_min{args.min_reads}_{lctag}.csv"
    active_x_csv = Path(args.active_x_csv).expanduser().resolve() if args.active_x_csv else default_active_x_csv
    haplotypes_csv = phase_dir / f"haplotypes_min{args.min_reads}_{lctag}.csv"

    if args.phased_vcf:
        requested_vcf = Path(args.phased_vcf).expanduser().resolve()
        if str(requested_vcf).endswith(".gz"):
            phased_vcf = requested_vcf
            phase_vcf_plain = Path(str(requested_vcf)[:-3])
        else:
            phase_vcf_plain = requested_vcf
            phased_vcf = Path(str(requested_vcf) + ".gz")
    else:
        phase_vcf_plain = phase_dir / f"haplotypes_min{args.min_reads}_{lctag}.vcf"
        phased_vcf = phase_vcf_plain.with_suffix(".vcf.gz")

    tagged_bam = Path(args.tagged_bam).expanduser().resolve() if args.tagged_bam else tagged_dir / f"{run_name}.tagged.bam"

    if not args.skip_matrices:
        if not args.original_gene_matrix_dir:
            raise ValueError("--original-gene-matrix-dir is required unless --skip-matrices is used")
        if not args.original_transcript_matrix_dir:
            raise ValueError("--original-transcript-matrix-dir is required unless --skip-matrices is used")
        original_gene_matrix_dir = Path(args.original_gene_matrix_dir).expanduser().resolve()
        original_transcript_matrix_dir = Path(args.original_transcript_matrix_dir).expanduser().resolve()
        ensure_dir(original_gene_matrix_dir, "Original gene matrix directory")
        ensure_dir(original_transcript_matrix_dir, "Original transcript matrix directory")
    else:
        original_gene_matrix_dir = None
        original_transcript_matrix_dir = None

    # 1. Count SNP reads.
    if args.skip_count:
        ensure_file(counts_tsv, "Existing counts TSV")
    else:
        if not maybe_skip("1. Count SNP reads", [counts_tsv], args.force):
            cmd = module_cmd(args.python, "count") + [
                "--bam", str(bam),
                "--vcf", str(vcf),
                "--out", str(counts_tsv),
                "--gtf", str(gtf),
                "--summary-out", str(count_summary),
                "--slop", str(args.slop),
                "--max-nearest-dist", str(args.max_nearest_dist),
                "--cb-tag", args.cb_tag,
                "--min-mapq", str(args.min_mapq),
                "--min-bq", str(args.min_bq),
                "--max-depth", str(args.max_depth),
                "--window", str(args.window),
                "--procs", str(args.procs),
                "--threads", str(args.threads),
            ]
            if args.region:
                cmd.extend(["--region", args.region])
            elif args.chrom:
                cmd.extend(["--chrom", args.chrom])
            add_opt(cmd, "--sample", args.count_sample)
            add_opt(cmd, "--read-evidence-out", args.read_evidence_out)
            add_flag(cmd, "--het-only", args.het_only)
            add_flag(cmd, "--allow-secondary", args.allow_secondary)
            add_flag(cmd, "--count-duplicates", args.count_duplicates)
            add_flag(cmd, "--skip-indels", args.skip_indels)
            add_flag(cmd, "--ignore-overlaps", args.ignore_overlaps)
            add_flag(cmd, "--drop-conflicts", args.drop_conflicts)
            add_flag(cmd, "--drop-multi-tsv-and-gtf", args.drop_multi_tsv_and_gtf)
            add_flag(cmd, "--debug", args.count_debug)

            wants_dropped = (args.drop_conflicts or args.drop_multi_tsv_and_gtf) and not args.no_dropped_read_outputs
            if wants_dropped:
                cmd.extend([
                    "--dropped-read-ids-out", str(dropped_read_evidence),
                    "--unique-dropped-read-ids-out", str(unique_dropped_read_ids),
                    "--dropped-read-summary-out", str(dropped_read_summary),
                ])

            run_step("1. Count SNP reads", cmd, commands_log, args.dry_run)

    # 2. Phase haplotypes.
    if args.skip_phase:
        ensure_file(active_x_csv, "Existing active-X CSV")
        ensure_file(phased_vcf, "Existing phased VCF")
    else:
        if not maybe_skip("2. Phase haplotypes", [haplotypes_csv, active_x_csv, phased_vcf], args.force):
            cmd = module_cmd(args.python, "phase") + [
                str(counts_tsv),
                str(phase_dir),
                "--min_reads", str(args.min_reads),
                "--lower_cutoff", str(args.lower_cutoff),
                "--partition_mode", args.partition_mode,
                "--vcf-out", str(phase_vcf_plain),
                "--sample-name", sample_name,
                "--stage2-mode", args.stage2_mode,
                "--stage2-min-evidence", str(args.stage2_min_evidence),
                "--stage2-tie-action", args.stage2_tie_action,
            ]
            add_opt(cmd, "--upper_cutoff", args.upper_cutoff)
            add_opt(cmd, "--edge_weight_mode", args.edge_weight_mode)
            add_opt(cmd, "--discordant_tsv", args.discordant_tsv)
            run_step("2. Phase haplotypes", cmd, commands_log, args.dry_run)

            if args.active_x_csv and not args.dry_run and default_active_x_csv.exists():
                active_x_csv.parent.mkdir(parents=True, exist_ok=True)
                shutil.copy2(default_active_x_csv, active_x_csv)

    # 3. Tag BAM.
    if args.skip_tag:
        ensure_file(tagged_bam, "Existing tagged BAM")
    else:
        if not maybe_skip("3. Tag BAM", [tagged_bam], args.force):
            cmd = module_cmd(args.python, "tag") + [
                str(bam),
                str(phased_vcf),
                str(active_x_csv),
                str(tagged_bam),
                "--mode", args.tag_mode,
                "--min-base-qual", str(args.tag_min_base_qual),
                "--min-score-diff", str(args.tag_min_score_diff),
                "--min-total-score", str(args.tag_min_total_score),
                "--min-count-diff", str(args.tag_min_count_diff),
                "--min-total-count", str(args.tag_min_total_count),
            ]
            run_step("3. Tag BAM", cmd, commands_log, args.dry_run)

    # 4. Split tagged BAM.
    expected_split_bams = [split_dir / f"{split_prefix}.{label}.bam" for label in ["X1", "X2", "Xa", "Xi", "low", "amb", "unknown"]]
    if args.skip_split:
        for pth in expected_split_bams:
            ensure_file(pth, "Existing split BAM")
    else:
        if not maybe_skip("4. Split tagged BAM", expected_split_bams, args.force):
            cmd = module_cmd(args.python, "split") + [
                str(tagged_bam),
                str(split_dir),
                "--prefix", split_prefix,
            ]
            run_step("4. Split tagged BAM", cmd, commands_log, args.dry_run)

    # 5. Build matrices.
    if not args.skip_matrices:
        expected_mats = [matrix_dir / "Xa.csv", matrix_dir / "Xi.csv"]
        if not maybe_skip("5. Build Xa/Xi matrices", expected_mats, args.force):
            cmd = module_cmd(args.python, "matrix") + [
                str(split_dir),
                str(original_gene_matrix_dir),
                str(original_transcript_matrix_dir),
                str(gtf),
                "--cb-tag", args.cb_tag,
                "--gn-tag", args.gn_tag,
                "--tr-tag", args.tr_tag,
                "--transcript-id-type", args.transcript_id_type,
                "--outdir", str(matrix_dir),
                "--bam-pattern", args.matrix_bam_pattern,
            ]
            if args.matrix_exclude_dropped_reads and unique_dropped_read_ids.exists():
                cmd.extend(["--exclude-read-ids", str(unique_dropped_read_ids)])
            add_flag(cmd, "--exclude-duplicates", args.matrix_exclude_duplicates)
            add_flag(cmd, "--keep-secondary", args.matrix_keep_secondary)
            add_flag(cmd, "--keep-supplementary", args.matrix_keep_supplementary)
            add_flag(cmd, "--skip-transcript", args.skip_transcript)
            run_step("5. Build Xa/Xi matrices", cmd, commands_log, args.dry_run)

    print("\nDone.")
    print(f"Commands log: {commands_log}")
    print("Main outputs:")
    print(f"  Counts TSV:      {counts_tsv}")
    print(f"  Active-X CSV:    {active_x_csv}")
    print(f"  Phased VCF:      {phased_vcf}")
    print(f"  Tagged BAM:      {tagged_bam}")
    print(f"  Split BAM dir:   {split_dir}")
    if not args.skip_matrices:
        print(f"  Matrix dir:      {matrix_dir}")

    return 0


def main(argv=None) -> int:
    args = parse_run_args(argv)
    return run_pipeline(args)


if __name__ == "__main__":
    raise SystemExit(main())
