
# scDaisychain
<img width="2126" height="2363" alt="fig1_less_verbose_outputs_network" src="https://github.com/user-attachments/assets/f2b68d4a-9533-40bf-b8c2-79ab610fb14a" />


scDaisychain is a toolkit for haplotype-resolved analysis of X chromosome inactivation single-cell long-read RNA sequencing data.

The software was developed and tested for Oxford Nanopore 10X single-cell datasets. Unlike short read 10X single cell sequencing where only the first 100 or 150bp of a transcript are typically captured, the longer read length of 10X nanopore sequencing allows for a greater number of SNPs to be used.

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
Which will install the python package depdencies `pandas`, `numpy`, `scipy`, `pysam`, `ipykernel`, `scanpy` and `python-igraph`.
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
or with curl:
```bash
curl -L \
  -o scDaisychain_mouse_chrX_test_data.tar.gz \
  https://github.com/daisy-kavanagh/scDaisychain/releases/download/v0.1.0/scDaisychain_mouse_chrX_test_data.tar.gz

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

## Downstream analysis
### Loading Matrices into scanpy
As scDaisychain produces Xa and Xi gene x cell matrices, they can easily be loaded as layers into your single cell analysis software of choice. We provide built in functions for loading them in scanpy. The haplotype sums from stage 1 of the phasing model file may also relevant to load for XCI skew analysis as it contains the active X call per cell.
```python
import scdaisychain as dc

matrix_dir = "/home/913/dk4874/scratch/gdata/scDaisychain_paper/scDaisychain/scDaisychain/test_data/gene_raw_feature_bc_matrix"
xa_csv = "/home/913/dk4874/scratch/gdata/scDaisychain_paper/mouse/output/scDaisychain_modes_multi_filtered_weighted_from_WGS/daisychain_split_bams_count/matrices/Xa.csv"
xi_csv = "/home/913/dk4874/scratch/gdata/scDaisychain_paper/mouse/output/scDaisychain_modes_multi_filtered_weighted_from_WGS/daisychain_split_bams_count/matrices/Xi.csv"
haplotype_sums_csv = "/home/913/dk4874/scratch/gdata/scDaisychain_paper/mouse/output/scDaisychain_modes_multi_filtered_weighted_from_WGS/haplotype_sums_df_min10_lc0.01.csv"

mouse = dc.load_10x_batches_with_optional_layers(
    matrices=[(matrix_dir, "mouse_BM")],
    sample_by_batch={
        "mouse_BM": "mouse_BM",
    },
    xa_csv_by_batch={
        "mouse_BM": xa_csv,
    },
    xi_csv_by_batch={
        "mouse_BM": xi_csv,
    },
    haplotype_sums_csv=haplotype_sums_csv,
    var_names="gene_symbols",
    make_dense_layers=False,   # safer for memory
)
```
This will produce an anndata object with the full gene expression layers loaded, an active_X observation for each cell stating the active X, and layers Xa, Xi. It also contains the layer Xai which is the sum of the Xa + Xi layers. Finally, it produces the layers 'Xa_Normalized', 'Xi_Normalized', 'Xai_Normalized', which are the respective X layers but divided by the total read count of the cell x 1000. 
Standard single cell processing such as filtering, clustering and cell type identification can then be proceeded with. 

After downstream clustering and cell type identification has been run, the per cell type per donor beta shrinkage can be applied. This can be run on the anndata object as follows:
```python
res = dc.run_xi_downstream_from_adata(
    adata,
    annotation="mouse_mouse_chrX_PAR_escape_annotation_for_xi_pipeline.tsv",
    outdir="xi_outputs",
    xi_layer="Xi",
    xa_layer="Xa",
    donor_col="sample",
    celltype_col="cluster_majority_celltype",
    gene_name_col="var_names",
)
```

Alternatively, if you have a h5ad file with the Xa/Xi layers loaded and cell type column it can be run from command line:
```bash
python xi_downstream_from_anndata.py \
  --h5ad adata.h5ad \
  --annotation mouse_mouse_chrX_PAR_escape_annotation_for_xi_pipeline.tsv \
  --outdir xi_outputs \
  --xi-layer Xi \
  --xa-layer Xa \
  --donor-col sample \
  --celltype-col cluster_majority_celltype \
  --gene-name-col var_names \
  --healthy-only false \
  --chrX-only true \
  --escape-mode recovery \
  --xi-cutoff 0.10 \
  --min-donors 1
```

