#!/bin/bash

#SBATCH -p general
#SBATCH -N 1
#SBATCH -n 2
#SBATCH -t 3-00:00:00
#SBATCH --mem=12G

# 检查是否提供了 SRA 文件列表
if [ "$#" -ne 1 ]; then
  echo "Usage: $0 <sra_list_file>"
  exit 1
fi

filelist=$1

# 定义工作目录
output_dir="output" #default output dir name
fq_dir="output/fastq"
bam_dir='output/bam'
snv='output/snv'
cfa_dir='output/cfa'
forup_dir='output/forup'
ano_dir='output/ano'
#snippy_snv='snippy_snv'

DIRS=("output/fastq" "output/bam" "output/cfa" "output/snv" "output/forup" "output/ano")
for dir in "${DIRS[@]}"; do
    if [ ! -d "$dir" ]; then
      mkdir -p "$dir"
    fi
done

# 检查样本是否已成功完成的函数
is_sample_done() {
  local sample=$1
  if [ -s "${snv}/${sample}.snp" ]; then
    return 0
  else
    return 1
  fi
}

# 处理每个 SRA 样本
cat $filelist |while read sample; do 
  # 检查样本是否已完成
  if is_sample_done "$sample"; then
    echo "Sample ${sample} is already completed. Skipping..."
    continue
  fi

  # 检查是否存在 fastq 文件
  if [ -f ${fq_dir}/${sample}_1.fastq.gz ] || [ -f ${fq_dir}/${sample}.fastq.gz ]; then
    echo "Fastq files for ${sample} exist. Skipping SRA extraction..."
  else
    echo "Fastq file for ${sample} doesn't exist."
    continue
  fi
  module purge
  # 进行后续处理
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
    cfa=${cfa_dir}/${sample}.cfa
    
    # 执行 paired-end 处理
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
    python ./src/remove_low_ebr.py $var > $ppe
    perl ./src/varscan_work_flow/1_format_trans.pl $ppe > $format
    perl ./src/varscan_work_flow/2_fix_extract.pl $format > $fix
    perl ./src/varscan_work_flow/3.1_mix_pileup_merge.pl $format $pileup > $forup
    cut -f2-4 $fix > $snp
    perl ./src/single_colony/1st_loci_recall_cns.pl $cns > $cfa
    #snp annotation
    perl ./src/mtbc_translate/1_MTBC_Annotation_mtbc_4411532.pl $snp > ${ano_dir}/${sample}.ano
    
    ##############call unfix mutation
    echo -e \"$sample\\c\" > kept/${sample}.cns.cov.dep
    sed 's/:/\t/g' $cns|awk '{if ($6 >= 3){n++;sum+=$6}} END {print "\t",n/4411532,"\t",sum/n}' >> kept/${sample}.cns.cov.dep
    #echo "NA" > kept/repeat_above20.list
    perl src/unfix_script/0_mix_extract_0.99.pl forup/${sample}.forup > kept/${sample}.mix
    perl src/unfix_script/1_forup_format.pl kept/${sample}.mix > kept/${sample}.mixfor
    perl src/unfix_script/2_info_mark.pl kept/${sample}.mixfor > kept/${sample}.mixmark
    perl src/unfix_script/3_redepin_filt.pl kept/repeat_above20.list kept/${sample}.cns.cov.dep kept/${sample}.mixmark
    rm kept/${sample}.cns kept/${sample}.cns.cov.dep kept/${sample}.mix kept/${sample}.mixfor kept/${sample}.mixmark kept/${sample}.mixmarkdisc
            
    #annotation
    awk '{print $9"\t"$10"\t"$11}' kept/${sample}.mixmarkkept > unfix_ano/${sample}.snp
    perl ./src/mtbc_translate/1_MTBC_Annotation_mtbc_4411532.pl unfix_ano/${sample}.snp > unfix_ano/${sample}_unfix.ano
    rm unfix_ano/${sample}.snp
    rm $sfq $fq1_tr $fq2_tr $fq3_tr $sai1 $sai2 $sai3 $samp $sams $bamp $bams $bamm $pileup $ppe $format $fix
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
    cfa=${cfa_dir}/${sample}.cfa
    
    # 执行 single-end 处理
    sickle se -f $sfq -t sanger -o $fq_tr
    bwa aln -t 2 -R 1 ./data/tb.ancestor.fasta $fq_tr > $sai
    bwa samse ./data/tb.ancestor.fasta $sai $fq_tr > $samf
    samtools view -bhSt ./data/tb.ancestor.fasta.fai $samf -o $bamf
    samtools sort $bamf -o $sortbam
    samtools index $sortbam
    samtools mpileup -q 30 -Q 20 -ABOf ./data/tb.ancestor.fasta $sortbam > $pileup
    java -jar ./src/VarScan.v2.3.9.jar mpileup2snp $pileup --min-coverage 3 --min-reads2 2 --min-avg-qual 20 --min-var-freq 0.01 --min-freq-for-hom 0.9 --p-value 99e-02 --strand-filter 0 > $var
    java -jar ./src/VarScan.v2.3.9.jar mpileup2cns $pileup --min-coverage 3 --min-avg-qual 20 --min-var-freq 0.75 --min-reads2 2 --strand-filter 0 > $cns
    python ./src/remove_low_ebr.py $var > $ppe
    perl ./src/varscan_work_flow/1_format_trans.pl $ppe > $format
    perl ./src/varscan_work_flow/2_fix_extract.pl $format > $fix
    perl ./src/varscan_work_flow/3.1_mix_pileup_merge.pl $format $pileup > $forup
    cut -f2-4 $fix > $snp
    #snp annotation
    perl ./src/mtbc_translate/1_MTBC_Annotation_mtbc_4411532.pl $snp > ${ano_dir}/${sample}.ano

    #perl ./src/single_colony/1st_loci_recall_cns.pl ./src/single_colony/all_4411532_loc.txt ./src/single_colony/depth.txt $cns > $cfa
    ###########call unfix mutation
    echo -e \"$sample\\c\" > kept/${sample}.cns.cov.dep
    sed 's/:/\t/g' $cns|awk '{if ($6 >= 3){n++;sum+=$6}} END {print "\t",n/4411532,"\t",sum/n}' >> kept/${sample}.cns.cov.dep
    #echo "NA" > kept/repeat_above20.list
    perl ./src/unfix_script/0_mix_extract_0.99.pl forup/${sample}.forup > kept/${sample}.mix
    perl ./src/unfix_script/1_forup_format.pl kept/${sample}.mix > kept/${sample}.mixfor
    perl ./src/unfix_script/2_info_mark.pl kept/${sample}.mixfor > kept/${sample}.mixmark
    perl ./src/unfix_script/3_redepin_filt.pl kept/repeat_above20.list kept/${sample}.cns.cov.dep kept/${sample}.mixmark
    rm kept/${sample}.cns kept/${sample}.cns.cov.dep kept/${sample}.mix kept/${sample}.mixfor kept/${sample}.mixmark kept/${sample}.mixmarkdisc
            
    #unfix annotation
    awk '{print $9"\t"$10"\t"$11}' kept/${sample}.mixmarkkept > unfix_ano/${sample}.snp
    perl ./src/mtbc_translate/1_MTBC_Annotation_mtbc_4411532.pl unfix_ano/${sample}.snp > unfix_ano/${sample}_unfix.ano
    rm unfix_ano/${sample}.snp
    rm $sfq $fq_tr $sai $samf $bamf $pileup $ppe $format $fix
    rm ${snv}/${sample}.cns.gz 
    gzip $cns
    echo "Processing complete for ${sample}."
  fi

  # 标记样本为已完成
  touch "${fq_dir}/${sample}_done"
done

echo "All samples processed."

