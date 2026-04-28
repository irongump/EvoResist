#!/bin/bash
#
# EvoResist Pipeline for candidate lists selection and optimization
# Selects optimal convergence event threshold or promoter length, followed
# by full-cohort evaluation and incremental-gain evaluation.
#
# Run from the script directory:
# 
#   bash 00-pipeline.sh
#

# ============================================================
# Path Configuration – edit these variables before running.
# ============================================================

DATA_DIR=./data          # input data folder (SNPs, indels, per-drug sample lists)
RESULTS_DIR=             # base results directory; per-drug outputs written here
SNP_DIR=                 # directory containing per-sample SNP annotation files
INDEL_DIR=               # directory containing per-sample indel annotation files

SNPs_100k_path=${DATA_DIR}/all_ann_convergent_flt.txt #this file should be the "{OUTDIR}/lineage_ann/all_ann_convergent_flt.txt"
indels_100k_path=${DATA_DIR}/all_indel_100k.txt

module load r
module load python

# ============================================================
# Generate 70/30 train/test ID splits per drug
# Input:  ./data/${drug_name}_sample_list.txt
# Output: ${RESULTS_DIR}/${drug_name}/id/train_70.txt
#         ${RESULTS_DIR}/${drug_name}/id/test_30.txt
# ============================================================

for drug_name in RIF INH EMB PZA LFX MFX BDQ AMK STM ETO KAN CAP LZD; do
    Rscript ./01-train_test_split.R \
        ${DATA_DIR}/${drug_name}_sample_list.txt \
        ${RESULTS_DIR}/${drug_name}/id/
done


# ============================================================
# Filter SNPs and indels from the 100k dataset to
# drug-relevant genomic regions
# ============================================================

filter_snp_allpos_inside() {
    local input_file="$1"
    local out_file="$2"
    local region_expr="$3"

    awk -F' ' -v REGION_EXPR="$region_expr" '
    $4 ~ /^[0-9]+(-[0-9]+)*$/ {
        n = split($4, a, "-")
        keep = 1
        for (i = 1; i <= n; i++) {
            pos = a[i]
            ok = 0
            if (REGION_EXPR == "RIF") {
                ok = (pos >= 758807 && pos <= 763325)
            } else if (REGION_EXPR == "INH") {
                ok = ((pos >= 1672440 && pos <= 1675011) || (pos >= 2153889 && pos <= 2157111))
            } else if (REGION_EXPR == "EMB") {
                ok = (pos >= 4245514 && pos <= 4249810)
            } else if (REGION_EXPR == "PZA") {
                ok = (pos >= 2288681 && pos <= 2290241)
            } else if (REGION_EXPR == "LFX_MFX") {
                ok = (pos >= 4240 && pos <= 9818)
            } else if (REGION_EXPR == "BDQ") {
                ok = ((pos >= 777990 && pos <= 779487) || (pos >= 1460045 && pos <= 1461290) || (pos >= 2858300 && pos <= 2861418))
            } else if (REGION_EXPR == "AMK") {
                ok = ((pos >= 1470846 && pos <= 1473382) || (pos >= 2715332 && pos <= 2716332))
            } else if (REGION_EXPR == "STM") {
                ok = ((pos >= 780560 && pos <= 781934) || (pos >= 4407528 && pos <= 4409202) || pos == 1472359 || pos == 1472362)
            } else if (REGION_EXPR == "ETO") {
                ok = ((pos >= 4326004 && pos <= 4328473) || (pos >= 1672440 && pos <= 1675011))
            } else if (REGION_EXPR == "KAN") {
                ok = ((pos >= 1470846 && pos <= 1473382) || (pos >= 2714124 && pos <= 2716332))
            } else if (REGION_EXPR == "CAP") {
                ok = ((pos >= 1470846 && pos <= 1473382) || (pos >= 1916940 && pos <= 1918746))
            } else if (REGION_EXPR == "LZD") {
                ok = ((pos >= 799809 && pos <= 801462) || (pos >= 1472658 && pos <= 1476795))
            }
            if (!ok) { keep = 0; break }
        }
        if (keep) print
    }
    ' "$input_file" > "$out_file"
}

### RIF
filter_snp_allpos_inside "${SNPs_100k_path}" ${RESULTS_DIR}/RIF/denovo_snp_2.txt "RIF"
awk -F'\t' '($2 >= 759807 && $2 <= 763325)' ${indels_100k_path} > ${RESULTS_DIR}/RIF/denovo_indel_2.txt

### INH
filter_snp_allpos_inside "${SNPs_100k_path}" ${RESULTS_DIR}/INH/denovo_snp_2.txt "INH"
awk -F'\t' '(($2 >= 1674202 && $2 <= 1675011) || ($2 >= 2153889 && $2 <= 2156111))' ${indels_100k_path} > ${RESULTS_DIR}/INH/denovo_indel_2.txt

### EMB
filter_snp_allpos_inside "${SNPs_100k_path}" ${RESULTS_DIR}/EMB/denovo_snp_2.txt "EMB"
awk -F'\t' '($2 >= 4246514 && $2 <= 4249810)' ${indels_100k_path} > ${RESULTS_DIR}/EMB/denovo_indel_2.txt

### PZA
filter_snp_allpos_inside "${SNPs_100k_path}" ${RESULTS_DIR}/PZA/denovo_snp_2.txt "PZA"
awk -F'\t' '($2 >= 2288681 && $2 <= 2289241)' ${indels_100k_path} > ${RESULTS_DIR}/PZA/denovo_indel_2.txt

### LFX
filter_snp_allpos_inside "${SNPs_100k_path}" ${RESULTS_DIR}/LFX/denovo_snp_2.txt "LFX_MFX"
awk -F'\t' '(($2 >= 5240 && $2 <= 7262) || ($2 >= 7302 && $2 <= 9818))' ${indels_100k_path} > ${RESULTS_DIR}/LFX/denovo_indel_2.txt

### MFX
filter_snp_allpos_inside "${SNPs_100k_path}" ${RESULTS_DIR}/MFX/denovo_snp_2.txt "LFX_MFX"
awk -F'\t' '(($2 >= 5240 && $2 <= 7262) || ($2 >= 7302 && $2 <= 9818))' ${indels_100k_path} > ${RESULTS_DIR}/MFX/denovo_indel_2.txt

### BDQ
filter_snp_allpos_inside "${SNPs_100k_path}" ${RESULTS_DIR}/BDQ/denovo_snp_2.txt "BDQ"
awk -F'\t' '(($2 >= 778990 && $2 <= 779487) || ($2 >= 1461045 && $2 <= 1461290) || ($2 >= 2859300 && $2 <= 2860418))' ${indels_100k_path} > ${RESULTS_DIR}/BDQ/denovo_indel_2.txt

### AMK
filter_snp_allpos_inside "${SNPs_100k_path}" ${RESULTS_DIR}/AMK/denovo_snp_2.txt "AMK"
awk -F'\t' '(($2 >= 1471846 && $2 <= 1473382) || ($2 >= 2715332 && $2 <= 2716332))' ${indels_100k_path} > ${RESULTS_DIR}/AMK/denovo_indel_2.txt

### STM
filter_snp_allpos_inside "${SNPs_100k_path}" ${RESULTS_DIR}/STM/denovo_snp_2.txt "STM"
awk -F'\t' '(($2 >= 781560 && $2 <= 781934) || ($2 >= 4407528 && $2 <= 4408202))' ${indels_100k_path} > ${RESULTS_DIR}/STM/denovo_indel_2.txt

### ETO
filter_snp_allpos_inside "${SNPs_100k_path}" ${RESULTS_DIR}/ETO/denovo_snp_2.txt "ETO"
awk -F'\t' '(($2 >= 4326004 && $2 <= 4327473) || ($2 >= 1674202 && $2 <= 1675011))' ${indels_100k_path} > ${RESULTS_DIR}/ETO/denovo_indel_2.txt

### KAN
filter_snp_allpos_inside "${SNPs_100k_path}" ${RESULTS_DIR}/KAN/denovo_snp_2.txt "KAN"
awk -F'\t' '(($2 >= 1471846 && $2 <= 1473382) || ($2 >= 2714124 && $2 <= 2716332))' ${indels_100k_path} > ${RESULTS_DIR}/KAN/denovo_indel_2.txt

### CAP
filter_snp_allpos_inside "${SNPs_100k_path}" ${RESULTS_DIR}/CAP/denovo_snp_2.txt "CAP"
awk -F'\t' '(($2 >= 1471846 && $2 <= 1473382) || ($2 >= 1917940 && $2 <= 1918746))' ${indels_100k_path} > ${RESULTS_DIR}/CAP/denovo_indel_2.txt

### LZD
filter_snp_allpos_inside "${SNPs_100k_path}" ${RESULTS_DIR}/LZD/denovo_snp_2.txt "LZD"
awk -F'\t' '(($2 >= 800809 && $2 <= 801462) || ($2 >= 1473658 && $2 <= 1476795))' ${indels_100k_path} > ${RESULTS_DIR}/LZD/denovo_indel_2.txt


# ============================================================
# Build annotated initial candidate variant list per drug
# Please check detailed information about input arguments within the R script.
# ============================================================

Rscript ./02-make_initial_list.R RIF rpoB 759807 763325 1 ${RESULTS_DIR}
Rscript ./02-make_initial_list.R INH inhA,katG 1674202,2153889 1675011,2156111 1,0 ${RESULTS_DIR}
Rscript ./02-make_initial_list.R EMB embB 4246514 4249810 1 ${RESULTS_DIR}
Rscript ./02-make_initial_list.R PZA pncA 2288681 2289241 0 ${RESULTS_DIR}
Rscript ./02-make_initial_list.R LFX gyrB,gyrA 5240,7302 7262,9818 1,1 ${RESULTS_DIR}
Rscript ./02-make_initial_list.R MFX gyrB,gyrA 5240,7302 7262,9818 1,1 ${RESULTS_DIR}
Rscript ./02-make_initial_list.R BDQ Rv0678,atpE,pepQ 778990,1461045,2859300 779487,1461290,2860418 1,1,0 ${RESULTS_DIR}
Rscript ./02-make_initial_list.R AMK rrs,eis 1471846,2714124 1473382,2715332 1,0 ${RESULTS_DIR}
Rscript ./02-make_initial_list.R STM gid,rpsL 4407528,781560 4408202,781934 0,1 ${RESULTS_DIR}
Rscript ./02-make_initial_list.R ETO inhA,ethA 1674202,4326004 1675011,4327473 1,0 ${RESULTS_DIR}
Rscript ./02-make_initial_list.R KAN rrs,eis 1471846,2714124 1473382,2715332 1,0 ${RESULTS_DIR}
Rscript ./02-make_initial_list.R CAP rrs,tlyA 1471846,1917940 1473382,1918746 1,1 ${RESULTS_DIR}
Rscript ./02-make_initial_list.R LZD rplC,rrl 800809,1473658 801462,1476795 1,1 ${RESULTS_DIR}




# ============================================================
# Threshold sweep (2-6) with fixed promoter = 500bp:
# generate list1 (convergence threshold + promoter filter applied)
# Please check detailed information about input arguments within the R script.
# ============================================================

promoter_region=500
for threshold in 2 3 4 5 6; do
    Rscript ./04-make_thres_prom_combination.R RIF rpoB 759807 763325 0 ${threshold} ${promoter_region} 1 ${RESULTS_DIR}
    Rscript ./04-make_thres_prom_combination.R INH inhA,katG 1674202,2153889 1675011,2156111 0,1 ${threshold} ${promoter_region} 1,0 ${RESULTS_DIR}
    Rscript ./04-make_thres_prom_combination.R EMB embB 4246514 4249810 0 ${threshold} ${promoter_region} 1 ${RESULTS_DIR}
    Rscript ./04-make_thres_prom_combination.R PZA pncA 2288681 2289241 1 ${threshold} ${promoter_region} 0 ${RESULTS_DIR}
    Rscript ./04-make_thres_prom_combination.R LFX gyrB,gyrA 5240,7302 7262,9818 0,0 ${threshold} ${promoter_region} 1,1 ${RESULTS_DIR}
    Rscript ./04-make_thres_prom_combination.R MFX gyrB,gyrA 5240,7302 7262,9818 0,0 ${threshold} ${promoter_region} 1,1 ${RESULTS_DIR}
    Rscript ./04-make_thres_prom_combination.R BDQ Rv0678,atpE,pepQ 778990,1461045,2859300 779487,1461290,2860418 1,0,1 ${threshold} ${promoter_region} 1,1,0 ${RESULTS_DIR}
    Rscript ./04-make_thres_prom_combination.R AMK rrs,eis 1471846,2714124 1473382,2715332 0,0 ${threshold} ${promoter_region} 1,0 ${RESULTS_DIR}
    Rscript ./04-make_thres_prom_combination.R STM gid,rpsL 4407528,781560 4408202,781934 1,0 ${threshold} ${promoter_region} 0,1 ${RESULTS_DIR}
    Rscript ./04-make_thres_prom_combination.R ETO inhA,ethA 1674202,4326004 1675011,4327473 0,1 ${threshold} ${promoter_region} 1,0 ${RESULTS_DIR}
    Rscript ./04-make_thres_prom_combination.R KAN rrs,eis 1471846,2714124 1473382,2715332 0,0 ${threshold} ${promoter_region} 1,0 ${RESULTS_DIR}
    Rscript ./04-make_thres_prom_combination.R CAP rrs,tlyA 1471846,1917940 1473382,1918746 0,1 ${threshold} ${promoter_region} 1,1 ${RESULTS_DIR}
    Rscript ./04-make_thres_prom_combination.R LZD rplC,rrl 800809,1473658 801462,1476795 0,0 ${threshold} ${promoter_region} 1,1 ${RESULTS_DIR}
done


# ============================================================
# Evaluate list1 on training set: leave-one-out (promoter = 500bp, threshold 2-6)
# ============================================================

promoter_region=500
for drug_name in RIF INH EMB PZA LFX MFX BDQ AMK STM ETO KAN CAP LZD; do
    for threshold in 2 3 4 5 6; do
        python3 ./03-leave_one_out.py \
            --variants_file ${RESULTS_DIR}/${drug_name}/Threshold_${threshold}_Promoter_${promoter_region}_list1.tsv \
            --id_list_file ${RESULTS_DIR}/${drug_name}/id/train_70.txt \
            --snp_dir ${SNP_DIR} \
            --indel_dir ${INDEL_DIR} \
            --output_dir ${RESULTS_DIR}/${drug_name}/Threshold_${threshold}_Promoter_${promoter_region}/train_list1/
    done
done


# ============================================================
# Apply leave-one-out variant removal criteria to generate list2
# (promoter = 500bp, threshold sweep 2-6)
# Please check detailed information about input arguments within the R script.
# ============================================================

promoter_region=500
for threshold in 2 3 4 5 6; do
    Rscript ./05-loo_evaluate.R RIF rpoB 759807 763325 0 ${threshold} ${promoter_region} 1 ${RESULTS_DIR}
    Rscript ./05-loo_evaluate.R INH inhA,katG 1674202,2153889 1675011,2156111 0,1 ${threshold} ${promoter_region} 1,0 ${RESULTS_DIR}
    Rscript ./05-loo_evaluate.R EMB embB 4246514 4249810 0 ${threshold} ${promoter_region} 1 ${RESULTS_DIR}
    Rscript ./05-loo_evaluate.R PZA pncA 2288681 2289241 1 ${threshold} ${promoter_region} 0 ${RESULTS_DIR}
    Rscript ./05-loo_evaluate.R LFX gyrB,gyrA 5240,7302 7262,9818 0,0 ${threshold} ${promoter_region} 1,1 ${RESULTS_DIR}
    Rscript ./05-loo_evaluate.R MFX gyrB,gyrA 5240,7302 7262,9818 0,0 ${threshold} ${promoter_region} 1,1 ${RESULTS_DIR}
    Rscript ./05-loo_evaluate.R BDQ Rv0678,atpE,pepQ 778990,1461045,2859300 779487,1461290,2860418 1,0,1 ${threshold} ${promoter_region} 1,1,0 ${RESULTS_DIR}
    Rscript ./05-loo_evaluate.R AMK rrs,eis 1471846,2714124 1473382,2715332 0,0 ${threshold} ${promoter_region} 1,0 ${RESULTS_DIR}
    Rscript ./05-loo_evaluate.R STM gid,rpsL 4407528,781560 4408202,781934 1,0 ${threshold} ${promoter_region} 0,1 ${RESULTS_DIR}
    Rscript ./05-loo_evaluate.R ETO inhA,ethA 1674202,4326004 1675011,4327473 0,1 ${threshold} ${promoter_region} 1,0 ${RESULTS_DIR}
    Rscript ./05-loo_evaluate.R KAN rrs,eis 1471846,2714124 1473382,2715332 0,0 ${threshold} ${promoter_region} 1,0 ${RESULTS_DIR}
    Rscript ./05-loo_evaluate.R CAP rrs,tlyA 1471846,1917940 1473382,1918746 0,1 ${threshold} ${promoter_region} 1,1 ${RESULTS_DIR}
    Rscript ./05-loo_evaluate.R LZD rplC,rrl 800809,1473658 801462,1476795 0,0 ${threshold} ${promoter_region} 1,1 ${RESULTS_DIR}
done


# ============================================================
# Evaluate list2 on train and test sets
# (promoter = 500bp, threshold sweep 2-6)
# ============================================================

promoter_region=500
for drug_name in RIF INH EMB PZA LFX MFX BDQ AMK STM ETO KAN CAP LZD; do
    for threshold in 2 3 4 5 6; do
        python3 ./03-leave_one_out.py \
            --variants_file ${RESULTS_DIR}/${drug_name}/Threshold_${threshold}_Promoter_${promoter_region}_list2.tsv \
            --id_list_file ${RESULTS_DIR}/${drug_name}/id/train_70.txt \
            --snp_dir ${SNP_DIR} \
            --indel_dir ${INDEL_DIR} \
            --output_dir ${RESULTS_DIR}/${drug_name}/Threshold_${threshold}_Promoter_${promoter_region}/train_list2/

        python3 ./03-leave_one_out.py \
            --variants_file ${RESULTS_DIR}/${drug_name}/Threshold_${threshold}_Promoter_${promoter_region}_list2.tsv \
            --id_list_file ${RESULTS_DIR}/${drug_name}/id/test_30.txt \
            --snp_dir ${SNP_DIR} \
            --indel_dir ${INDEL_DIR} \
            --output_dir ${RESULTS_DIR}/${drug_name}/Threshold_${threshold}_Promoter_${promoter_region}/test_list2/
    done
done


# ============================================================
# Fix best threshold (select from test set) per drug, sweep promoter lengths (100-1000bp):
# generate list1
# Best thresholds: RIF=6, INH=3, EMB=5, PZA=2, LFX=6, MFX=6,
#                  BDQ=3, AMK=6, STM=5, ETO=4, KAN=6, CAP=6, LZD=5
# ============================================================

for promoter_region in 100 200 300 400 600 700 800 900 1000; do
    Rscript ./04-make_thres_prom_combination.R RIF rpoB 759807 763325 0 6 ${promoter_region} 1 ${RESULTS_DIR}
    Rscript ./04-make_thres_prom_combination.R INH inhA,katG 1674202,2153889 1675011,2156111 0,1 3 ${promoter_region} 1,0 ${RESULTS_DIR}
    Rscript ./04-make_thres_prom_combination.R EMB embB 4246514 4249810 0 5 ${promoter_region} 1 ${RESULTS_DIR}
    Rscript ./04-make_thres_prom_combination.R PZA pncA 2288681 2289241 1 2 ${promoter_region} 0 ${RESULTS_DIR}
    Rscript ./04-make_thres_prom_combination.R LFX gyrB,gyrA 5240,7302 7262,9818 0,0 6 ${promoter_region} 1,1 ${RESULTS_DIR}
    Rscript ./04-make_thres_prom_combination.R MFX gyrB,gyrA 5240,7302 7262,9818 0,0 6 ${promoter_region} 1,1 ${RESULTS_DIR}
    Rscript ./04-make_thres_prom_combination.R BDQ Rv0678,atpE,pepQ 778990,1461045,2859300 779487,1461290,2860418 1,0,1 3 ${promoter_region} 1,1,0 ${RESULTS_DIR}
    Rscript ./04-make_thres_prom_combination.R AMK rrs,eis 1471846,2714124 1473382,2715332 0,0 6 ${promoter_region} 1,0 ${RESULTS_DIR}
    Rscript ./04-make_thres_prom_combination.R STM gid,rpsL 4407528,781560 4408202,781934 1,0 5 ${promoter_region} 0,1 ${RESULTS_DIR}
    Rscript ./04-make_thres_prom_combination.R ETO inhA,ethA 1674202,4326004 1675011,4327473 0,1 4 ${promoter_region} 1,0 ${RESULTS_DIR}
    Rscript ./04-make_thres_prom_combination.R KAN rrs,eis 1471846,2714124 1473382,2715332 0,0 6 ${promoter_region} 1,0 ${RESULTS_DIR}
    Rscript ./04-make_thres_prom_combination.R CAP rrs,tlyA 1471846,1917940 1473382,1918746 0,1 6 ${promoter_region} 1,1 ${RESULTS_DIR}
    Rscript ./04-make_thres_prom_combination.R LZD rplC,rrl 800809,1473658 801462,1476795 0,0 5 ${promoter_region} 1,1 ${RESULTS_DIR}
done


# ============================================================
# Evaluate list1 on training set (fixed threshold, promoter sweep)
# ============================================================

for promoter_region in 100 200 300 400 600 700 800 900 1000; do

    threshold=6
    for drug_name in AMK CAP KAN LFX MFX RIF; do
        python3 ./03-leave_one_out.py \
            --variants_file ${RESULTS_DIR}/${drug_name}/Threshold_${threshold}_Promoter_${promoter_region}_list1.tsv \
            --id_list_file ${RESULTS_DIR}/${drug_name}/id/train_70.txt \
            --snp_dir ${SNP_DIR} \
            --indel_dir ${INDEL_DIR} \
            --output_dir ${RESULTS_DIR}/${drug_name}/Threshold_${threshold}_Promoter_${promoter_region}/train_list1/
    done

    threshold=5
    for drug_name in EMB STM LZD; do
        python3 ./03-leave_one_out.py \
            --variants_file ${RESULTS_DIR}/${drug_name}/Threshold_${threshold}_Promoter_${promoter_region}_list1.tsv \
            --id_list_file ${RESULTS_DIR}/${drug_name}/id/train_70.txt \
            --snp_dir ${SNP_DIR} \
            --indel_dir ${INDEL_DIR} \
            --output_dir ${RESULTS_DIR}/${drug_name}/Threshold_${threshold}_Promoter_${promoter_region}/train_list1/
    done

    threshold=4
    for drug_name in ETO; do
        python3 ./03-leave_one_out.py \
            --variants_file ${RESULTS_DIR}/${drug_name}/Threshold_${threshold}_Promoter_${promoter_region}_list1.tsv \
            --id_list_file ${RESULTS_DIR}/${drug_name}/id/train_70.txt \
            --snp_dir ${SNP_DIR} \
            --indel_dir ${INDEL_DIR} \
            --output_dir ${RESULTS_DIR}/${drug_name}/Threshold_${threshold}_Promoter_${promoter_region}/train_list1/
    done

    threshold=3
    for drug_name in BDQ INH; do
        python3 ./03-leave_one_out.py \
            --variants_file ${RESULTS_DIR}/${drug_name}/Threshold_${threshold}_Promoter_${promoter_region}_list1.tsv \
            --id_list_file ${RESULTS_DIR}/${drug_name}/id/train_70.txt \
            --snp_dir ${SNP_DIR} \
            --indel_dir ${INDEL_DIR} \
            --output_dir ${RESULTS_DIR}/${drug_name}/Threshold_${threshold}_Promoter_${promoter_region}/train_list1/
    done

    threshold=2
    for drug_name in PZA; do
        python3 ./03-leave_one_out.py \
            --variants_file ${RESULTS_DIR}/${drug_name}/Threshold_${threshold}_Promoter_${promoter_region}_list1.tsv \
            --id_list_file ${RESULTS_DIR}/${drug_name}/id/train_70.txt \
            --snp_dir ${SNP_DIR} \
            --indel_dir ${INDEL_DIR} \
            --output_dir ${RESULTS_DIR}/${drug_name}/Threshold_${threshold}_Promoter_${promoter_region}/train_list1/
    done

done


# ============================================================
# Apply variant removal criteria to generate list2 (promoter sweep)
# ============================================================

for promoter_region in 100 200 300 400 600 700 800 900 1000; do
    Rscript ./05-loo_evaluate.R RIF rpoB 759807 763325 0 6 ${promoter_region} 1 ${RESULTS_DIR}
    Rscript ./05-loo_evaluate.R INH inhA,katG 1674202,2153889 1675011,2156111 0,1 3 ${promoter_region} 1,0 ${RESULTS_DIR}
    Rscript ./05-loo_evaluate.R EMB embB 4246514 4249810 0 5 ${promoter_region} 1 ${RESULTS_DIR}
    Rscript ./05-loo_evaluate.R PZA pncA 2288681 2289241 1 2 ${promoter_region} 0 ${RESULTS_DIR}
    Rscript ./05-loo_evaluate.R LFX gyrB,gyrA 5240,7302 7262,9818 0,0 6 ${promoter_region} 1,1 ${RESULTS_DIR}
    Rscript ./05-loo_evaluate.R MFX gyrB,gyrA 5240,7302 7262,9818 0,0 6 ${promoter_region} 1,1 ${RESULTS_DIR}
    Rscript ./05-loo_evaluate.R BDQ Rv0678,atpE,pepQ 778990,1461045,2859300 779487,1461290,2860418 1,0,1 3 ${promoter_region} 1,1,0 ${RESULTS_DIR}
    Rscript ./05-loo_evaluate.R AMK rrs,eis 1471846,2714124 1473382,2715332 0,0 6 ${promoter_region} 1,0 ${RESULTS_DIR}
    Rscript ./05-loo_evaluate.R STM gid,rpsL 4407528,781560 4408202,781934 1,0 5 ${promoter_region} 0,1 ${RESULTS_DIR}
    Rscript ./05-loo_evaluate.R ETO inhA,ethA 1674202,4326004 1675011,4327473 0,1 4 ${promoter_region} 1,0 ${RESULTS_DIR}
    Rscript ./05-loo_evaluate.R KAN rrs,eis 1471846,2714124 1473382,2715332 0,0 6 ${promoter_region} 1,0 ${RESULTS_DIR}
    Rscript ./05-loo_evaluate.R CAP rrs,tlyA 1471846,1917940 1473382,1918746 0,1 6 ${promoter_region} 1,1 ${RESULTS_DIR}
    Rscript ./05-loo_evaluate.R LZD rplC,rrl 800809,1473658 801462,1476795 0,0 5 ${promoter_region} 1,1 ${RESULTS_DIR}
done


# ============================================================
# Evaluate list2 on train and test sets (fixed threshold, promoter sweep)
# ============================================================

for promoter_region in 100 200 300 400 600 700 800 900 1000; do

    threshold=6
    for drug_name in AMK CAP KAN LFX MFX RIF; do
        python3 ./03-leave_one_out.py \
            --variants_file ${RESULTS_DIR}/${drug_name}/Threshold_${threshold}_Promoter_${promoter_region}_list2.tsv \
            --id_list_file ${RESULTS_DIR}/${drug_name}/id/train_70.txt \
            --snp_dir ${SNP_DIR} \
            --indel_dir ${INDEL_DIR} \
            --output_dir ${RESULTS_DIR}/${drug_name}/Threshold_${threshold}_Promoter_${promoter_region}/train_list2/

        python3 ./03-leave_one_out.py \
            --variants_file ${RESULTS_DIR}/${drug_name}/Threshold_${threshold}_Promoter_${promoter_region}_list2.tsv \
            --id_list_file ${RESULTS_DIR}/${drug_name}/id/test_30.txt \
            --snp_dir ${SNP_DIR} \
            --indel_dir ${INDEL_DIR} \
            --output_dir ${RESULTS_DIR}/${drug_name}/Threshold_${threshold}_Promoter_${promoter_region}/test_list2/
    done

    threshold=5
    for drug_name in EMB STM LZD; do
        python3 ./03-leave_one_out.py \
            --variants_file ${RESULTS_DIR}/${drug_name}/Threshold_${threshold}_Promoter_${promoter_region}_list2.tsv \
            --id_list_file ${RESULTS_DIR}/${drug_name}/id/train_70.txt \
            --snp_dir ${SNP_DIR} \
            --indel_dir ${INDEL_DIR} \
            --output_dir ${RESULTS_DIR}/${drug_name}/Threshold_${threshold}_Promoter_${promoter_region}/train_list2/

        python3 ./03-leave_one_out.py \
            --variants_file ${RESULTS_DIR}/${drug_name}/Threshold_${threshold}_Promoter_${promoter_region}_list2.tsv \
            --id_list_file ${RESULTS_DIR}/${drug_name}/id/test_30.txt \
            --snp_dir ${SNP_DIR} \
            --indel_dir ${INDEL_DIR} \
            --output_dir ${RESULTS_DIR}/${drug_name}/Threshold_${threshold}_Promoter_${promoter_region}/test_list2/
    done

    threshold=4
    for drug_name in ETO; do
        python3 ./03-leave_one_out.py \
            --variants_file ${RESULTS_DIR}/${drug_name}/Threshold_${threshold}_Promoter_${promoter_region}_list2.tsv \
            --id_list_file ${RESULTS_DIR}/${drug_name}/id/train_70.txt \
            --snp_dir ${SNP_DIR} \
            --indel_dir ${INDEL_DIR} \
            --output_dir ${RESULTS_DIR}/${drug_name}/Threshold_${threshold}_Promoter_${promoter_region}/train_list2/

        python3 ./03-leave_one_out.py \
            --variants_file ${RESULTS_DIR}/${drug_name}/Threshold_${threshold}_Promoter_${promoter_region}_list2.tsv \
            --id_list_file ${RESULTS_DIR}/${drug_name}/id/test_30.txt \
            --snp_dir ${SNP_DIR} \
            --indel_dir ${INDEL_DIR} \
            --output_dir ${RESULTS_DIR}/${drug_name}/Threshold_${threshold}_Promoter_${promoter_region}/test_list2/
    done

    threshold=3
    for drug_name in BDQ INH; do
        python3 ./03-leave_one_out.py \
            --variants_file ${RESULTS_DIR}/${drug_name}/Threshold_${threshold}_Promoter_${promoter_region}_list2.tsv \
            --id_list_file ${RESULTS_DIR}/${drug_name}/id/train_70.txt \
            --snp_dir ${SNP_DIR} \
            --indel_dir ${INDEL_DIR} \
            --output_dir ${RESULTS_DIR}/${drug_name}/Threshold_${threshold}_Promoter_${promoter_region}/train_list2/

        python3 ./03-leave_one_out.py \
            --variants_file ${RESULTS_DIR}/${drug_name}/Threshold_${threshold}_Promoter_${promoter_region}_list2.tsv \
            --id_list_file ${RESULTS_DIR}/${drug_name}/id/test_30.txt \
            --snp_dir ${SNP_DIR} \
            --indel_dir ${INDEL_DIR} \
            --output_dir ${RESULTS_DIR}/${drug_name}/Threshold_${threshold}_Promoter_${promoter_region}/test_list2/
    done

    threshold=2
    for drug_name in PZA; do
        python3 ./03-leave_one_out.py \
            --variants_file ${RESULTS_DIR}/${drug_name}/Threshold_${threshold}_Promoter_${promoter_region}_list2.tsv \
            --id_list_file ${RESULTS_DIR}/${drug_name}/id/train_70.txt \
            --snp_dir ${SNP_DIR} \
            --indel_dir ${INDEL_DIR} \
            --output_dir ${RESULTS_DIR}/${drug_name}/Threshold_${threshold}_Promoter_${promoter_region}/train_list2/

        python3 ./03-leave_one_out.py \
            --variants_file ${RESULTS_DIR}/${drug_name}/Threshold_${threshold}_Promoter_${promoter_region}_list2.tsv \
            --id_list_file ${RESULTS_DIR}/${drug_name}/id/test_30.txt \
            --snp_dir ${SNP_DIR} \
            --indel_dir ${INDEL_DIR} \
            --output_dir ${RESULTS_DIR}/${drug_name}/Threshold_${threshold}_Promoter_${promoter_region}/test_list2/
    done

done


# ============================================================
# Final evaluation on the full per-drug sample list
# Evaluates EvoResist final list, WHO G1+G2, and WHO G1 lists.
# ============================================================

for drug_name in RIF INH EMB PZA LFX MFX BDQ AMK STM ETO KAN CAP LZD; do
    python3 ./03-leave_one_out.py \
        --variants_file ${RESULTS_DIR}/${drug_name}/EvoResist_Final_list.tsv \
        --id_list_file ${DATA_DIR}/${drug_name}_sample_list.txt \
        --snp_dir ${SNP_DIR} \
        --indel_dir ${INDEL_DIR} \
        --output_dir ${RESULTS_DIR}/${drug_name}/Final_Evaluate/EvoResist/

    python3 ./03-leave_one_out.py \
        --variants_file ${RESULTS_DIR}/${drug_name}/WHO_G1_2_withcolnames.txt \
        --id_list_file ${DATA_DIR}/${drug_name}_sample_list.txt \
        --snp_dir ${SNP_DIR} \
        --indel_dir ${INDEL_DIR} \
        --output_dir ${RESULTS_DIR}/${drug_name}/Final_Evaluate/G1_2/

    python3 ./03-leave_one_out.py \
        --variants_file ${RESULTS_DIR}/${drug_name}/WHO_G1_withcolnames.txt \
        --id_list_file ${DATA_DIR}/${drug_name}_sample_list.txt \
        --snp_dir ${SNP_DIR} \
        --indel_dir ${INDEL_DIR} \
        --output_dir ${RESULTS_DIR}/${drug_name}/Final_Evaluate/G1/
done


# ============================================================
# Evaluate incremental gain of newly prioritized mutations
# relative to the WHO G1+G2 baseline
# ============================================================

for drug_name in RIF INH EMB PZA LFX MFX BDQ LZD AMK STM ETO KAN CAP; do
    python3 ./06-gain_evaluation.py \
        --list1_variants_file ${RESULTS_DIR}/${drug_name}/WHO_G1_2_withcolnames.txt \
        --list2_variants_file ${RESULTS_DIR}/${drug_name}/EvoResist_Final_list.tsv \
        --id_list_file ${DATA_DIR}/${drug_name}_sample_list.txt \
        --snp_dir ${SNP_DIR} \
        --indel_dir ${INDEL_DIR} \
        --output_file ${RESULTS_DIR}/${drug_name}/gain_EvoResist.txt
done
