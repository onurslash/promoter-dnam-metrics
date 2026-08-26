#!/usr/bin/env bash
#
# 01_download_and_align.sh
#
# Raw data to per-CpG tables for the CD3+ analysis.
#
#   SRA download -> Trim Galore -> BSMAP (hg19) -> CAMDA.py -> per-CpG tables
#   SRA download -> STAR -> RSEM -> transcript FPKM
#
# Produces the three inputs consumed by 02_analysis.qmd:
#   data/Methylation.tsv        per-CpG methylation ratio  (BSMAP methratio schema)
#   data/CAMDA.tsv              per-CpG concurrence ratio  (same schema, C_count = concurrence)
#   data/CD3.isoforms.results   RSEM transcript quantification
#
# Reference files fetched here as well:
#   ref/hg19.fa, ref/gencode.v19.annotation.gtf, ref/refGene.txt, ref/cgi.bed
#
# Requirements: ~200 GB free disk, ~32 GB RAM, several hours of CPU time.
# CAMDA.py requires samtools 0.1.19 specifically; a separate conda environment
# is created for it because the flag it uses was removed in samtools 1.x.
#
# Usage:  bash 01_download_and_align.sh
# Env:    THREADS (default 8), ROOT (default $PWD)

set -euo pipefail

ROOT="${ROOT:-$PWD}"
THREADS="${THREADS:-8}"
WGBS_SRR="SRR1104838"     # GSM1186660, CD3+ T cells, 37-year-old male donor
RNA_SRR="SRR980468"       # GSM1220574, matched RNA-seq from the same donor

mkdir -p "$ROOT"/{data,ref,fastq,align,tmp}
cd "$ROOT"
export TMPDIR="$ROOT/tmp"

# ---------------------------------------------------------------- environments
if ! conda env list | grep -q '^cd3 '; then
  conda create -y -n cd3 -c conda-forge -c bioconda \
    sra-tools pigz fastqc trim-galore bsmap star=2.6.0c rsem \
    python=2.7 pysam r-base r-data.table
fi
if ! conda env list | grep -q '^cd3-samtools019 '; then
  conda create -y -n cd3-samtools019 -c bioconda samtools=0.1.19
fi
# shellcheck disable=SC1091
eval "$(conda shell.bash hook)"
conda activate cd3
SAMTOOLS019="$(conda run -n cd3-samtools019 which samtools)"

[[ -d CAMDA ]] || git clone https://github.com/JiejunShi/CAMDA.git

# ---------------------------------------------------------------- references
if [[ ! -s ref/hg19.fa ]]; then
  wget -q -O ref/hg19.fa.gz https://hgdownload.soe.ucsc.edu/goldenPath/hg19/bigZips/hg19.fa.gz
  gunzip ref/hg19.fa.gz
fi
if [[ ! -s ref/gencode.v19.annotation.gtf ]]; then
  wget -q -O ref/g.gtf.gz \
    https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_19/gencode.v19.annotation.gtf.gz
  gunzip -c ref/g.gtf.gz > ref/gencode.v19.annotation.gtf && rm ref/g.gtf.gz
fi
if [[ ! -s ref/refGene.txt ]]; then
  wget -q -O ref/refGene.txt.gz https://hgdownload.soe.ucsc.edu/goldenPath/hg19/database/refGene.txt.gz
  gunzip ref/refGene.txt.gz
fi
if [[ ! -s ref/cgi.bed ]]; then
  wget -q -O ref/cpgIslandExt.txt.gz \
    https://hgdownload.soe.ucsc.edu/goldenPath/hg19/database/cpgIslandExt.txt.gz
  gunzip -c ref/cpgIslandExt.txt.gz | awk 'BEGIN{OFS="\t"}{print $2,$3,$4}' > ref/cgi.bed
  rm ref/cpgIslandExt.txt.gz
fi

# ---------------------------------------------------------------- fastq
fetch () {
  local srr=$1
  [[ -s fastq/${srr}_1.fastq.gz ]] && return
  prefetch --max-size u "$srr"
  fasterq-dump --split-files -e "$THREADS" -t "$TMPDIR" -O fastq "$srr"
  pigz -p "$THREADS" fastq/"${srr}"_*.fastq
}
fetch "$WGBS_SRR"
fetch "$RNA_SRR"

# ---------------------------------------------------------------- WGBS
if [[ ! -s align/wgbs.bam ]]; then
  if [[ -s fastq/${WGBS_SRR}_2.fastq.gz ]]; then
    trim_galore --paired --gzip --cores 4 -o fastq \
      "fastq/${WGBS_SRR}_1.fastq.gz" "fastq/${WGBS_SRR}_2.fastq.gz"
    bsmap -a "fastq/${WGBS_SRR}_1_val_1.fq.gz" -b "fastq/${WGBS_SRR}_2_val_2.fq.gz" \
          -d ref/hg19.fa -o align/wgbs.bam -p "$THREADS"
  else
    trim_galore --gzip --cores 4 -o fastq "fastq/${WGBS_SRR}_1.fastq.gz"
    bsmap -a "fastq/${WGBS_SRR}_1_trimmed.fq.gz" \
          -d ref/hg19.fa -o align/wgbs.bam -p "$THREADS"
  fi
fi

# CAMDA.py emits the concurrence ratio and the conventional methylation ratio
# in a common per-CpG format; -x CG restricts output to CpG context.
python CAMDA/scripts/CAMDA.py CAMDA \
  align/wgbs.bam ref/hg19.fa -o data/CD3 -w data/CD3 -s "$SAMTOOLS019" -x CG

mv -f data/CD3_CpG_CAMDA.tsv     data/CAMDA.tsv
mv -f data/CD3_CpG_MethRatio.tsv data/Methylation.tsv

# ---------------------------------------------------------------- RNA-seq
if [[ ! -d ref/star ]]; then
  mkdir -p ref/star
  STAR --runMode genomeGenerate --genomeDir ref/star \
       --genomeFastaFiles ref/hg19.fa --sjdbGTFfile ref/gencode.v19.annotation.gtf \
       --sjdbOverhang 100 --runThreadN "$THREADS"
fi
if [[ ! -d ref/rsem ]]; then
  mkdir -p ref/rsem
  rsem-prepare-reference --gtf ref/gencode.v19.annotation.gtf ref/hg19.fa ref/rsem/hg19
fi

RNA2=""; [[ -s fastq/${RNA_SRR}_2.fastq.gz ]] && RNA2="fastq/${RNA_SRR}_2.fastq.gz"
STAR --genomeDir ref/star --readFilesIn "fastq/${RNA_SRR}_1.fastq.gz" $RNA2 \
     --readFilesCommand zcat --quantMode TranscriptomeSAM \
     --outSAMtype BAM Unsorted --outFileNamePrefix align/rna. --runThreadN "$THREADS"

PAIRED=""; [[ -n "$RNA2" ]] && PAIRED="--paired-end"
rsem-calculate-expression $PAIRED --bam --no-bam-output -p "$THREADS" \
  align/rna.Aligned.toTranscriptome.out.bam ref/rsem/hg19 data/CD3

mv -f data/CD3.isoforms.results data/CD3.isoforms.results

echo
echo "Done. Inputs for 02_analysis.qmd:"
ls -lh data/Methylation.tsv data/CAMDA.tsv data/CD3.isoforms.results
