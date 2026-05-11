#!/usr/bin/env bash
set -euo pipefail

###############################################################################
# EvoResist INDEL calling + QC
#
# This script performs:
#   1) Basic BAM filtering (remove unmapped; MAPQ>=1)
#   2) Add read groups (Picard)
#   3) GATK3 local realignment around indels (RealignerTargetCreator + IndelRealigner)
#   4) GATK3 UnifiedGenotyper indel calling
#   5) Indel normalization (bcftools norm; left-align + split multiallelic)
#   6) Add AF/AC/AN tags (bcftools +fill-tags)
#   7) QC filtering:
#        - keep INDELs only
#        - indel length <= 50 bp
#        - depth >= 10 (FORMAT/DP)
#        - allele fraction >= 0.75 (INFO/AF)
#
# Outputs per sample:
#   <sample>.INDEL.raw.vcf
#   <sample>.INDEL.norm.vcf.gz
#   <sample>.INDEL.norm.tags.vcf.gz
#   <sample>.INDEL.pass.vcf.gz  (+ .tbi index)
#
# Requirements:
#   - samtools
#   - java (JDK8 compatible with GATK3)
#   - picard.jar
#   - GATK3 GenomeAnalysisTK.jar
#   - bcftools (with +fill-tags)
#
# Notes:
#   - This is a "minimal modification" of the original pipeline: we append
#     normalization and QC steps after UnifiedGenotyper.
#   - Bidirectional (strand) support filtering is NOT enforced here because
#     UnifiedGenotyper VCFs may not contain strand-count tags (e.g., DP4).
###############################################################################

############################################
# User-configurable paths
############################################
REF="data/tb.ancestor.fasta"
PICARD_JAR="scripts/picard.jar" #install by users
GATK_JAR="scripts/GenomeAnalysisTK.jar" #install by users

############################################
# QC thresholds (match manuscript/rebuttal)
############################################
MAX_INDEL_LEN=50
MIN_DP=10
MIN_AF=0.75

############################################
# Input: directory + file pattern
# Usage:
#   bash indel_call_qc.sh /path/to/bams "*.sort.bam"
############################################
BAM_DIR="${1:-.}" # should be "{OUTDIR}/bam"
BAM_GLOB="${2:-*.sort.bam}"

# Optional: load modules (uncomment if needed)
# module load jdk/1.8.0_172-fasrc01
# module load samtools
# module load bcftools

cd "$BAM_DIR"

echo "[INFO] Working directory: $(pwd)"
echo "[INFO] Reference: $REF"
echo "[INFO] BAM pattern: $BAM_GLOB"
echo "[INFO] QC: indel_len<=${MAX_INDEL_LEN}, DP>=${MIN_DP}, AF>=${MIN_AF}"

# Sanity checks
command -v samtools >/dev/null 2>&1 || { echo "[ERROR] samtools not found"; exit 1; }
command -v bcftools >/dev/null 2>&1 || { echo "[ERROR] bcftools not found"; exit 1; }
[[ -s "$REF" ]] || { echo "[ERROR] Reference FASTA not found: $REF"; exit 1; }
[[ -s "$PICARD_JAR" ]] || { echo "[ERROR] Picard jar not found: $PICARD_JAR"; exit 1; }
[[ -s "$GATK_JAR" ]] || { echo "[ERROR] GATK jar not found: $GATK_JAR"; exit 1; }

# Loop over BAMs
shopt -s nullglob
bams=( $BAM_GLOB )
if [[ ${#bams[@]} -eq 0 ]]; then
  echo "[ERROR] No BAMs matched pattern: $BAM_GLOB"
  exit 1
fi

for bam in "${bams[@]}"; do
  base="$(basename "$bam")"
  prefix="${base%.sort.bam}"

  echo "============================================================"
  echo "[INFO] Sample: $prefix"
  echo "[INFO] Input BAM: $bam"
  echo "============================================================"

  # Intermediate files (same naming style as original pipeline)
  bam_filt="${prefix}.filt.bam"
  bam_rg="${prefix}.RG.bam"
  intervals="${prefix}.intervals"
  bam_rl="${prefix}.RL.bam"

  # VCFs
  vcf_raw="${prefix}.INDEL.raw.vcf"
  vcf_norm="${prefix}.INDEL.norm.vcf.gz"
  vcf_tags="${prefix}.INDEL.norm.tags.vcf.gz"
  vcf_pass="${prefix}.INDEL.pass.vcf.gz"

  ###########################################################################
  # 1) Filter BAM: remove unmapped reads (-F 4) and MAPQ=0 reads (-q 1)
  ###########################################################################
  samtools view -bF 4 -q 1 "$bam" > "$bam_filt"
  samtools index "$bam_filt"

  ###########################################################################
  # 2) Add or replace read groups (required by GATK tools)
  #    NOTE: RG fields are placeholders as in original script.
  ###########################################################################
  java -jar "$PICARD_JAR" AddOrReplaceReadGroups \
    I="$bam_filt" O="$bam_rg" \
    RGID=4 RGLB=lib1 RGPL=illumina RGPU=unit1 RGSM=20
  samtools index "$bam_rg"

  ###########################################################################
  # 3) GATK3: identify realignment targets around indels
  ###########################################################################
  java -Xmx2g -jar "$GATK_JAR" \
    -I "$bam_rg" -R "$REF" \
    -T RealignerTargetCreator -maxInterval 20000 \
    -o "$intervals"

  ###########################################################################
  # 4) GATK3: local realignment around indels
  ###########################################################################
  java -Xmx4g -jar "$GATK_JAR" \
    -I "$bam_rg" -R "$REF" \
    -T IndelRealigner -targetIntervals "$intervals" \
    -o "$bam_rl"
  samtools index "$bam_rl"

  ###########################################################################
  # 5) GATK3: call INDELs (UnifiedGenotyper)
  ###########################################################################
  java -Xmx4g -jar "$GATK_JAR" \
    -T UnifiedGenotyper -nt 1 \
    -R "$REF" -I "$bam_rl" \
    -glm INDEL -o "$vcf_raw"

  ###########################################################################
  # 6) Normalize INDEL representation (left-align; split multiallelic)
  ###########################################################################
  bcftools norm -f "$REF" -m -both "$vcf_raw" -Oz -o "$vcf_norm"
  bcftools index -t "$vcf_norm"

  ###########################################################################
  # 7) Add tags AF/AC/AN (computed from genotype fields)
  #    Note: DP is NOT supported by +fill-tags (expected behavior).
  ###########################################################################
  bcftools +fill-tags "$vcf_norm" -Oz -o "$vcf_tags" -- -t AF,AC,AN
  bcftools index -t "$vcf_tags"

  ###########################################################################
  # 8) QC filter:
  #    - keep only INDELs
  #    - indel length <= 50 bp
  #    - FORMAT/DP >= 10  (per-sample depth; avoids INFO/DP vs FORMAT/DP ambiguity)
  #    - INFO/AF >= 0.75  (fixed/high-confidence variants)
  ###########################################################################
  bcftools view -i "TYPE='indel' && abs(strlen(REF)-strlen(ALT))<=${MAX_INDEL_LEN}" "$vcf_tags" -Ou \
    | bcftools filter -i "FORMAT/DP>=${MIN_DP} && INFO/AF>=${MIN_AF}" -Oz -o "$vcf_pass"
  bcftools index -t "$vcf_pass"

  ###########################################################################
  # 9) Cleanup intermediates (keep key outputs)
  ###########################################################################
  rm -f "$bam_filt" "$bam_filt.bai" \
        "$bam_rg" "$bam_rg.bai" \
        "$intervals" \
        "$bam_rl" "$bam_rl.bai"

  echo "[INFO] Outputs:"
  echo "  - $vcf_raw"
  echo "  - $vcf_norm (+.tbi)"
  echo "  - $vcf_tags (+.tbi)"
  echo "  - $vcf_pass (+.tbi)"
done

echo "[INFO] Done."