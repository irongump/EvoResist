#!/bin/bash

#SBATCH -p general
#SBATCH -N 1
#SBATCH -n 2
#SBATCH -t 10-00:00:00
#SBATCH --mem=12G

# Check if a SRA file list is provided
if [ "$#" -ne 1 ]; then
  echo "Usage: $0 <sra_list_file>"
  exit 1
fi

filelist=$1

# Define working directories
output_dir="output" #default output dir name
fq_dir="output/fastq"
SRA="./data/sra"
bam_dir='output/bam'
snv='output/snv'
cfa_dir='output/cfa'
forup_dir='output/forup'
ano_dir='output/ano'

DIRS=("output/fastq" "output/bam" "output/cfa" "output/snv" "output/forup")
for dir in "${DIRS[@]}"; do
    if [ ! -d "$dir" ]; then
      mkdir -p "$dir"
    fi
done

# Function to check if a sample has already been processed
is_sample_done() {
  local sample=$1
  if [ -s "${snv}/${sample}.snp" ]; then
    return 0
  else
    return 1
  fi
}

# Process each SRA sample
cat $filelist |while read sample; do 
  # Check if sample is already completed
  if is_sample_done "$sample"; then
    echo "Sample ${sample} is already completed. Skipping..."
    continue
  fi

  module add sratoolkit/3.0.7
  # Check if FASTQ files already exist
  if [ -f ${fq_dir}/${sample}_1.fastq.gz ] || [ -f ${fq_dir}/${sample}.fastq.gz ]; then
    echo "Fastq files for ${sample} exist. Skipping SRA extraction..."
  else
    # Extract FASTQ from SRA files
    echo "Extracting fastq files for ${sample} from SRA..."
    if [ -f "${SRA}/${sample}/${sample}.sra" ]; then
      fastq-dump -split-3 --gzip -O ${fq_dir} ${SRA}/${sample}/${sample}.sra
    elif [ -f "${SRA}/${sample}/${sample}.sralite" ]; then
      fastq-dump -split-3 --gzip -O ${fq_dir} ${SRA}/${sample}/${sample}.sralite
    else
      echo "SRA file for ${sample} not found. Skipping..."
      continue
    fi
  fi
  module purge
  # Proceed with downstream processing
  fq1=${fq_dir}/${sample}_1.fastq.gz
  fq2=${fq_dir}/${sample}_2.fastq.gz
  sfq=${fq_dir}/${sample}.fastq.gz

  if [ -f "$fq2" ]; then
    echo "Processing paired-end reads for ${sample}..."
    fq1_tr=${fq_dir}/${sample}_tr_1.fq
    fq2_tr=${fq_dir}/${sample}_tr_2.fq
    fq3_tr=${fq_dir}/${sample}_tr_S.fq
    sai1=${fq_dir}/${sample}_R1.sai
    sai2=${fq_dir}/${sample}_R2.sai
    sai3=${fq_dir}/${sample}_S.sai
    samp=${bam_dir}/${sample}.paired.sam
    sams=${bam_dir}/${sample}.single.sam
    bamp=${bam_dir}/${sample}.paired.bam
    bams=${bam_dir}/${sample}.single.bam
    bamm=${bam_dir}/${sample}.merge.bam
    sortbam=${bam_dir}/${sample}.sort.bam
    pileup=${snv}/${sample}.pileup
    var=${snv}/${sample}.varscan
    cns=${snv}/${sample}.cns
    ppe=${snv}/${sample}.ppe
    format=${snv}/${sample}.for
    forup=${forup_dir}/${sample}.forup
    fix=${snv}/${sample}.fix
    snp=${snv}/${sample}.snp
    sample_cfa=${cfa_dir}/${sample}.cfa
    
    # Paired-end processing pipeline
    sickle pe -l 35 -f $fq1 -r $fq2 -t sanger -o $fq1_tr -p $fq2_tr -s $fq3_tr
    bwa aln -R 1 -t 2 ./data/tb.ancestor.fasta $fq1_tr > $sai1
    bwa aln -R 1 -t 2 ./data/tb.ancestor.fasta $fq2_tr > $sai2
    bwa aln -R 1 ./data/tb.ancestor.fasta $fq3_tr > $sai3
    bwa sampe -a 1000 -n 1 -N 1 ./data/tb.ancestor.fasta $sai1 $sai2 $fq1_tr $fq2_tr > $samp
    bwa samse -n 1 ./data/tb.ancestor.fasta $sai3 $fq3_tr > $sams
    samtools view -bhSt ./data/tb.ancestor.fasta.fai $samp -o $bamp
    samtools view -bhSt ./data/tb.ancestor.fasta.fai $sams -o $bams
    samtools merge $bamm $bamp $bams
    samtools sort $bamm -o $sortbam
    samtools index $sortbam
    samtools mpileup -q 30 -Q 20 -ABOf ./data/tb.ancestor.fasta $sortbam > $pileup
    java -jar ./src/VarScan.v2.3.9.jar mpileup2snp $pileup --min-coverage 3 --min-reads2 2 --min-avg-qual 20 --min-var-freq 0.01 --min-freq-for-hom 0.9 --p-value 99e-02 --strand-filter 0 > $var
    java -jar ./src/VarScan.v2.3.9.jar mpileup2cns $pileup --min-coverage 3 --min-avg-qual 20 --min-var-freq 0.75 --min-reads2 2 --strand-filter 0 > $cns
    python scripts/remove_low_ebr.py ./data/RLC_lowmapK50E4_H37Rv_pos.txt $var > $ppe
    perl scripts/varscan_work_flow/1_format_trans.pl $ppe > $format
    perl scripts/varscan_work_flow/2_fix_extract.pl $format > $fix
    perl scripts/varscan_work_flow/3.1_mix_pileup_merge.pl $format $pileup > $forup
    cut -f2-4 $fix > $snp
    perl scripts/single_colony/1st_loci_recall_cns.pl $cns > $sample_cfa
    rm $fq1 $fq2 $sfq $fq1_tr $fq2_tr $fq3_tr $sai1 $sai2 $sai3 $samp $sams $bamp $bams $bamm $pileup $ppe $format $fix
    rm ${snv}/${sample}.cns.gz 
    gzip $cns
    echo "Processing complete for ${sample}."
  else
    echo "Processing single-end reads for ${sample}..."
    fq_tr=${fq_dir}/${sample}_tr.fq
    sai=${fq_dir}/${sample}.sai
    samf=${bam_dir}/${sample}.sam
    bamf=${bam_dir}/${sample}.bam
    sortbam=${bam_dir}/${sample}.sort.bam
    pileup=${snv}/${sample}.pileup
    var=${snv}/${sample}.varscan
    cns=${snv}/${sample}.cns
    ppe=${snv}/${sample}.ppe
    format=${snv}/${sample}.for
    fix=${snv}/${sample}.fix
    snp=${snv}/${sample}.snp
    forup=${forup_dir}/${sample}.forup
    sample_cfa=${cfa_dir}/${sample}.cfa
    
    # Single-end processing pipeline
    sickle se -f $sfq -t sanger -o $fq_tr
    bwa aln -t 2 -R 1 ./data/tb.ancestor.fasta $fq_tr > $sai
    bwa samse ./data/tb.ancestor.fasta $sai $fq_tr > $samf
    samtools view -bhSt ./data/tb.ancestor.fasta.fai $samf -o $bamf
    samtools sort $bamf -o $sortbam
    samtools index $sortbam
    samtools mpileup -q 30 -Q 20 -ABOf ./data/tb.ancestor.fasta $sortbam > $pileup
    java -jar ./src/VarScan.v2.3.9.jar mpileup2snp $pileup --min-coverage 3 --min-reads2 2 --min-avg-qual 20 --min-var-freq 0.01 --min-freq-for-hom 0.9 --p-value 99e-02 --strand-filter 0 > $var
    java -jar ./src/VarScan.v2.3.9.jar mpileup2cns $pileup --min-coverage 3 --min-avg-qual 20 --min-var-freq 0.75 --min-reads2 2 --strand-filter 0 > $cns
    python scripts/remove_low_ebr.py ./data/RLC_lowmapK50E4_H37Rv_pos.txt $var > $ppe
    perl scripts/varscan_work_flow/1_format_trans.pl $ppe > $format
    perl scripts/varscan_work_flow/2_fix_extract.pl $format > $fix
    perl scripts/varscan_work_flow/3.1_mix_pileup_merge.pl $format $pileup > $forup
    cut -f2-4 $fix > $snp
    perl scripts/single_colony/1st_loci_recall_cns.pl $cns > $sample_cfa
    rm $sfq $fq_tr $sai $samf $bamf $pileup $ppe $format $fix
    rm ${snv}/${sample}.cns.gz 
    gzip $cns
    echo "Processing complete for ${sample}."
  fi

  # Mark sample as completed
  touch "${fq_dir}/${sample}_done"
done

echo "All samples processed."
