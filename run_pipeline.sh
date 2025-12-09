#!/usr/bin/env bash
set -euo pipefail

# Usage / argument handling

if [ "$#" -lt 4 ]; then
  echo "Usage: $0 REF.fasta R1.fastq R2.fastq OUT_PREFIX"
  exit 1
fi

REF=$1    # reference genome FASTA
R1=$2     # reaad 1 fastq (Illumina)
R2=$3     # read 2 fastq (Illumina)
PREFIX=$4 # prefix for all outputs, e.g. results/part2/EcoliK12

mkdir -p "$(dirname "$PREFIX")"

# 1. Map reads and sort BAM (minimap2+samtools)

echo "[RUN] Mapping reads with minimap2 and sorting BAM"

if [ "$R2" = "-" ]; then
  echo "[RUN] Single-end mode (no R2)"
  minimap2 -a -x sr "$REF" "$R1" \
    | samtools view -h -F 0x900 - \
    | samtools sort -O bam -o "${PREFIX}.sorted.bam"
else
  echo "[RUN] Paired-end mode"
  minimap2 -a -x sr "$REF" "$R1" "$R2" \
    | samtools view -h -F 0x900 - \
    | samtools sort -O bam -o "${PREFIX}.sorted.bam"
fi

samtools index "${PREFIX}.sorted.bam"


# 2. Variant caller 1 - bcftools

echo "[RUN] Calling variants with bcftools"

bcftools mpileup -Ou -f "$REF" "${PREFIX}.sorted.bam" \
  | bcftools call -mv -Oz -o "${PREFIX}.bcftools.vcf.gz"
  
bcftools index "${PREFIX}.bcftools.vcf.gz"


# 3. Variant caller 2 - Snippy

echo "[RUN] Calling variants with snippy"

SNIPPY_DIR="${PREFIX}_snippy_out"
SNIPPY_VCF="${SNIPPY_DIR}/snps.vcf"

if [ "$R2" = "-" ]; then
  echo "[RUN] Snippy single-end mode"
  snippy --cpus 4 \
         --outdir "$SNIPPY_DIR" \
         --ref "$REF" \
         --se "$R1" \
         --force
else
  echo "[RUN] Snippy paired-end mode"
  snippy --cpus 4 \
         --outdir "$SNIPPY_DIR" \
         --ref "$REF" \
         --R1 "$R1" \
         --R2 "$R2" \
         --force
fi

if [ ! -f "$SNIPPY_VCF" ]; then
    echo "[ERROR] Snippy did not produce $SNIPPY_VCF - cannot continue."
    exit 1
fi

# 4. Normalise VCFs (split multiallelic, left-align)

echo "[RUN] Normalising VCFs with bcftools norm"

bcftools norm -f "$REF" -m-any "${PREFIX}.bcftools.vcf.gz" -Oz -o "${PREFIX}.bcftools.norm.vcf.gz"
bcftools norm -f "$REF" -m-any "$SNIPPY_VCF" -Oz -o "${PREFIX}.snippy.norm.vcf.gz"

bcftools index "${PREFIX}.bcftools.norm.vcf.gz"
bcftools index "${PREFIX}.snippy.norm.vcf.gz"

# 5. Add INFO tag to say which caller produced which variant

echo "[RUN] Tagging variants with CALLER=bcftools/snippy"

bcftools annotate -x INFO/CALLER -I +'%CHROM:%POS:%REF:%ALT,bcftools' \
  "${PREFIX}.bcftools.norm.vcf.gz" \
  -Oz -o "${PREFIX}.bcftools.tagged.vcf.gz"

bcftools annotate -x INFO/CALLER -I +'%CHROM:%POS:%REF:%ALT,snippy' \
  "${PREFIX}.snippy.norm.vcf.gz" \
  -Oz -o "${PREFIX}.snippy.tagged.vcf.gz"

bcftools index "${PREFIX}.bcftools.tagged.vcf.gz"
bcftools index "${PREFIX}.snippy.tagged.vcf.gz"

# 6. Combine results of both callers into 1 VCF

echo "[RUN] Creating combined VCF (bcftools + snippy, no genotype columns)"

# Drop all sample/genotype columns (-G) so sample-name conflicts disappear
bcftools view -G "${PREFIX}.bcftools.tagged.vcf.gz" -Oz -o "${PREFIX}.bcftools.tagged.nosample.vcf.gz"
bcftools view -G "${PREFIX}.snippy.tagged.vcf.gz"   -Oz -o "${PREFIX}.snippy.tagged.nosample.vcf.gz"

bcftools index "${PREFIX}.bcftools.tagged.nosample.vcf.gz"
bcftools index "${PREFIX}.snippy.tagged.nosample.vcf.gz"

# Now concat the two no-sample VCFs
bcftools concat -a \
  "${PREFIX}.bcftools.tagged.nosample.vcf.gz" \
  "${PREFIX}.snippy.tagged.nosample.vcf.gz" \
  -Oz -o "${PREFIX}.combined.vcf.gz"

bcftools index "${PREFIX}.combined.vcf.gz"

echo "[DONE] Combined VCF: ${PREFIX}.combined.vcf.gz"


