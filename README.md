# scDaisychain

Installable package for the XCI workflow:

1. Count per-cell REF/ALT SNP reads from BAM + VCF.
2. Phase X-linked SNPs with scDaisychain graph phasing.
3. Tag reads with H1/H2 and Xa/Xi assignments.
4. Split the tagged BAM into X1/X2/Xa/Xi/low/amb/unknown BAMs.
5. Build gene/transcript x cell matrices from the split BAMs.

## Install

Recommended on HPC:

```bash
conda env create -f environment.yml
conda activate scDaisychain
```

Or, from an existing environment:

```bash
pip install -e .
```

Python package dependencies are `pandas`, `numpy`, `scipy`, `pysam`, and `python-igraph`.
The phasing step also calls external `bgzip` and `bcftools index`, so install `bcftools`/`htslib` through conda or your module system.

## Commands

The package installs one full-pipeline command and one command per step:

```bash
scDaisychain --help
scDaisychain run --help
xci-count-snps --help
xci-phase-x --help
xci-tag-bam --help
xci-split-bam --help
xci-make-matrices --help
```

## Full pipeline example

```bash
scDaisychain run \
  --bam input.bam \
  --vcf variants.vcf.gz \
  --gtf genes.gtf.gz \
  --outdir scdaisychain_run \
  --original-gene-matrix-dir /path/to/gene_matrix \
  --original-transcript-matrix-dir /path/to/transcript_matrix \
  --chrom chrX \
  --het-only \
  --min-reads 10 \
  --lower-cutoff 0.01 \
  --partition-mode expression \
  --tag-mode count
```

The output layout is:

```text
scdaisychain_run/
├── 01_counts/
├── 02_phase/
├── 03_tagged_bam/
├── 04_split_bams/
├── 05_matrices/
└── logs/
    └── commands.sh
```

Use `--dry-run` to print the exact commands without running them:

```bash
scDaisychain run ... --dry-run
```

## Step-by-step commands

```bash
xci-count-snps --help
xci-phase-x --help
xci-tag-bam --help
xci-split-bam --help
xci-make-matrices --help
```

Each module is import-safe and follows the pattern:

```python
def parse_args(argv=None):
    ...

def main(argv=None):
    ...

if __name__ == "__main__":
    raise SystemExit(main())
```
