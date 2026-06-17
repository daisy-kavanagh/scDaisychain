# scDaisychain

scDaisychain is a toolkit for haplotype-resolved analysis of X chromosome inactivation single-cell long-read RNA sequencing data.

The software was developed and tested for Oxford Nanopore 10X single-cell datasets. Unlike short read 10X single cell sequencing where only the first 100bp of a transcript are typically captured, the longer read length of 10X nanopore sequencing allows for a greater number of SNPs to be used.

Features

- Allelic quantification of SNPs at single cell resolution
- Phasing of X chromosome SNPs based on single cell expression
- Active X call per cell for XCI skew analysis
- Matrices of expression from the active and inactive X
- Matrices of expression from the two parental X haplotypes
- Downstream analysis tools for loading into scanpy, identifying XCI escapees and XCI skew at a celltype level.

Requirements

- 10X nanopore single cell data
- A vcf of X chromosome variants present in the sample (recommended to be generated from DNA sequencing)

## Install

First clone the github repo
```bash
git clone https://github.com/daisy-kavanagh/scDaisychain.git
cd scDaisychain
```

Easiest way to install is by creating a new mamba / conda environment which will include all dependencies:

```bash

mamba env create -f environment.yml
mamba activate scDaisychain
```

Or you can install an existing environment:

```bash
pip install -e .
```
Which will install the python package depdencies `pandas`, `numpy`, `scipy`, `pysam`, `scanpy` and `python-igraph`.
The phasing step also calls external `bgzip` and `bcftools`, so install `bcftools`/`htslib` through mamba/conda or your module system.


## Preprocessing
Data Preprocessing The software requires a BAM file of single cell nanopore reads, with tags for gene name (GN), cell barcode (CB) and UMI (UB), and the raw gene expression matrices. This can be achieved for example with the [epi2me wf-singlecell pipeline](https://github.com/epi2me-labs/wf-single-cell).

The resulting BAM must be deduplicated for example with [UMI tools](https://github.com/CGATOxford/UMI-tools). For efficiency, you can filter to just chrX with samtools first, as only a chrX bam is needed for the downstream steps.

An unphased vcf of variants is required also. For example for short read DNA sequencing with [GATK haplotype caller](https://gatk.broadinstitute.org/hc/en-us/articles/360037225632-HaplotypeCaller) or for Nanopore WGS [Clair3](https://github.com/HKU-BAL/Clair3).


## Example usage

Example mouse chrX C57/Bl6 X CAST/EiJ test data are available from the GitHub Releases page and can be downloaded as follows.

```bash
wget https://github.com/daisy-kavanagh/scDaisychain/releases/download/v0.1.0/scDaisychain_mouse_chrX_test_data.tar.gz
tar -xzf scDaisychain_mouse_chrX_test_data.tar.gz
```
Which will produce a test_data folder with the following files:
```text
test_data/
├── CASTxBL6_F1.0p1.vcf.gz
├── CASTxBL6_F1.0p1.vcf.gz.tbi
├── deduped.chrX.bam
├── deduped.chrX.bam.bai
├── gene_raw_feature_bc_matrix
│   ├── barcodes.tsv.gz
│   ├── features.tsv.gz
│   └── matrix.mtx.gz
├── genes.gtf
└── transcript_raw_feature_bc_matrix
    ├── barcodes.tsv.gz
    ├── features.tsv.gz
    └── matrix.mtx.gz
```

The full scDaisychain pipeline can be run with one command. The run parameters can be viewed with:
```bash
scDaisychain run --help
```

The following is an example of how to run the pipeline using the provided mouse data:
```bash
scDaisychain run \
  --bam test_data/deduped.chrX.bam \
  --vcf test_data/CASTxBL6_F1.0p1.vcf.gz \
  --gtf test_data/genes.gtf \
  --outdir scdaisychain_run \
  --original-gene-matrix-dir test_data/gene_raw_feature_bc_matrix \
  --original-transcript-matrix-dir test_data/transcript_raw_feature_bc_matrix \
  --min-reads 10 \
  --lower-cutoff 0.01 \
  --partition-mode weighted \
  --tag-mode count \
  --drop-conflicts \
  --drop-multi-tsv-and-gtf
```

The output layout is:

```text
scdaisychain_run/
├── 01_counts/ #Per cell allele counts of the variants in the VCF
├── 02_phase/ #Outputs of the phasing algorithm. Includes a phased chrX vcf, and a haplotype_sums_df.csv, that contains the active X call for each cell.
├── 03_tagged_bam/ #The input BAM with reads tagged as Xa or Xi (active/inactive) and X1/X2 (parental1/2)
├── 04_split_bams/ #BAM files split by the haplotype tags from the tagged bam stage
├── 05_matrices/ #Gene expression matrices from the split BAMs which can be used for downstream analysis of XCI escape
└── logs/
    └── commands.sh
```

Individual modules of the pipeline can also be run separately. Their help information can be viewed as follows:
```bash
scDaisychain-count-snps --help
scDaisychain-phase-x --help
scDaisychain-tag-bam --help
scDaisychain-split-bam --help
scDaisychain-make-matrices --help
```



