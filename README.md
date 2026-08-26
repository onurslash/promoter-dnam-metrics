# Promoter DNA Methylation Metrics

Reproducible code for the manuscript **“Promoter DNA Methylation Metrics and Gene Expression in Human CD3+ T Cells: Sensitivity to Analytical Choices.”**

This repository reproduces the analysis comparing regional mean methylation, the CAMDA concurrence ratio, and single-CpG β-entropy in matched WGBS and RNA-seq from human CD3+ T cells. The main aim is to evaluate how promoter methylation–expression associations change with analytical choices such as annotation source, transcript- versus gene-level quantification, coverage threshold, regional aggregation, treatment of zero-expression transcripts, promoter position, and CpG-island definition.

## Repository structure

```text
.
├── 01_download_and_align.sh   # download, alignment, methylation/CAMDA calling, RNA-seq quantification
├── 02_analysis.qmd            # complete statistical analysis, tables, bootstrap analyses and figures
├── data/                      # generated analysis inputs (not intended for version control)
├── ref/                       # downloaded reference files
├── fastq/                     # downloaded FASTQ files
├── align/                     # intermediate alignment files
├── figures/                   # figures produced by 02_analysis.qmd
└── tmp/                       # temporary files
```

## Public data

The analysis uses matched data from a single human CD3+ T-cell donor:

- **WGBS:** GEO `GSM1186660`, SRA `SRR1104838`
- **RNA-seq:** GEO `GSM1220574`, SRA `SRR980468`
- **Genome build:** hg19
- **Primary annotation:** GENCODE v19
- **Sensitivity annotation:** UCSC `refGene`
- **CpG islands:** UCSC hg19 `cpgIslandExt`

All required reference files are downloaded automatically by `01_download_and_align.sh`.

## Workflow

The workflow has two stages.

### 1. Download and process the sequencing data

```bash
bash 01_download_and_align.sh
```

The script performs:

```text
WGBS:
SRA -> Trim Galore -> BSMAP -> CAMDA.py -> per-CpG methylation/CAMDA tables

RNA-seq:
SRA -> STAR -> RSEM -> transcript-level FPKM
```

It produces the three files used by the analysis:

```text
data/Methylation.tsv
data/CAMDA.tsv
data/CD3.isoforms.results
```

The number of threads and working directory can be changed with environment variables:

```bash
THREADS=16 ROOT=/path/to/project bash 01_download_and_align.sh
```

### 2. Reproduce the statistical analysis

Render the Quarto document from the repository root:

```bash
quarto render 02_analysis.qmd
```

`02_analysis.qmd` computes the reported results directly from the processed data and generates the manuscript tables and figures. The analysis includes:

- promoter windows from −1,000 to +500 bp around the strand-corrected TSS;
- fifteen non-overlapping 100 bp bins;
- regional mean methylation, CAMDA concurrence ratio, and single-CpG β-entropy;
- pooled and CpG-weighted whole-promoter aggregation;
- fixed-promoter positional analysis;
- paired promoter bootstrap and chromosome-block bootstrap;
- the 24 estimable configurations of the analytical sensitivity grid;
- common-support marginal effects of analytical choices;
- expression-dependent and expression-independent representative-transcript rules;
- CpG-island stratification and alternative island definitions;
- CpG-density and sequencing-coverage sensitivity analyses;
- aggregation-scale analysis of the two regional entropy formulations.

A fixed random seed is used for bootstrap analyses.

## Main figure outputs

The Quarto workflow writes figures to `figures/` in both PNG and PDF format:

```text
Figure1_bin_profile
Figure2_delta_bootstrap
Figure3_full_grid
Figure4_cgi_strata
Figure5_aggregation_scale
```

## Computational requirements

The complete workflow is intended for a Linux environment with Conda available.

Approximate requirements for the raw-data workflow:

- **Disk:** ~200 GB free space
- **RAM:** ~32 GB recommended
- **CPU time:** several hours, depending on hardware

The downstream Quarto analysis uses large per-CpG tables (approximately 2.3 GB each), reaches roughly 16 GB peak memory, and typically takes about 1–2 hours on a multicore workstation.

### Important legacy dependency

The reference CAMDA implementation requires **samtools 0.1.19**. `01_download_and_align.sh` therefore creates a separate Conda environment for this version. The main processing environment contains the remaining tools, including SRA Toolkit, Trim Galore, BSMAP, STAR, RSEM, Python 2.7, R, and `data.table`.

Quarto must also be installed and available on `PATH` to render `02_analysis.qmd`.

## Reproducibility notes

- Methylation and CAMDA values are evaluated on matched CpG sets within each comparison.
- CpGs failing the specified coverage filter are removed before regional aggregation.
- Zero-expression transcripts are treated explicitly as an analytical factor rather than silently removed in every analysis.
- Transcript-level RNA-seq quantification is based on GENCODE v19; RefSeq is therefore used only for gene-level annotation comparisons.
- The sensitivity design contains **24 estimable configurations**, rather than a nominal fully crossed 32-cell design, because RefSeq transcript-level abundance is not directly estimable from the GENCODE-based RNA-seq quantification.
- Generated FASTQ, BAM, reference and large per-CpG files are not required to be stored in the GitHub repository because they can be regenerated from the public accessions above.

## Software provenance

The CAMDA implementation is obtained by the processing script from the original repository:

- Shi et al., *Nature Communications* (2021), **The concurrence of DNA methylation and demethylation is associated with transcription regulation**. DOI: `10.1038/s41467-021-25521-7`

Other software versions and the R session used for the statistical analysis are recorded by the scripts and by `sessionInfo()` in the rendered Quarto output.

## License

This repository is released under the MIT License. See `LICENSE` for details.

## Author

**R. Onur Öztornacı**  
Department of Biostatistics and Medical Informatics  
Faculty of Medicine, Gaziantep University, Türkiye
