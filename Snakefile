"""
EvoResist Snakemake Pipeline (v2.0)
====================================
Evolution-Guided Prioritization of Drug Resistance Mutations Enhances
Molecular Prediction of Tuberculosis Drug Susceptibility.

Pipeline Steps (Steps 1–7: convergent evolution analysis; Steps 8–14: DR mutation selection):
  1. snp_calling      - SNP calling from FASTQ or SRA files (per sample)
  2. indel_calling    - INDEL calling from BAM files (per sample)
  3. build_tree       - Phylogenetic tree building with IQ-TREE (per lineage/sublineage)
  4. branch_mutations - Branch mutations extraction (within lineage/sublineage)
  5. ancestor_mutations - Ancestor mutations extraction (prior to lineage/sublineage diversification)
  6. merge_annotations / stat_convergent / filter_convergent - Count convergent mutations by codon
  7. simulation       - GTR+Gamma simulation of mutations under a null distribution
  8. dr_train_test_split  - Stratify 70/30 train-test for sensitivity analysis
  9. dr_filter_variants   - Include convergent SNPs or indels from DR genes and promoter regions
  10. dr_initial_list      - Build annotated initial candidate variant list per drug
  11. dr_make_list1         - Apply threshold × promoter-length combination to generate list1
  12. dr_loo_evaluate       - Leave-one-out evaluation of a variant list on train or test split
  13. dr_make_list2         - Apply LOO filtering criteria to generate refined list2
  14. (manual) pending      - Final evaluation and incremental gain analysis

Key configuration options (config/config.yaml):
  outdir      - Base output directory (default: "output").
                Recommended: single-level relative path so that shell scripts
                using ../../ back-references resolve correctly.
  stop_at     - Stop after a specific step; valid values:
                  step1/snp_calling, step2/indel_calling,
                  step3/build_tree, step4/branch_mutations,
                  step5/ancestor_mutations, step6/merge_annotations,
                  step7/simulation/all,
                  step8/dr_prep, step9/dr_selection
  input_type  - Input data format: auto (default), fastq, or sra.
                When set to "fastq", fastq_dir must also be set.
  fastq_dir   - Directory with pre-existing per-sample FASTQ files.
                Required when input_type is "fastq"; ignored otherwise.
  lineage     - Process only the specified lineage or list of lineages
                (e.g. "Lineage1.1" or ["Lineage1.1", "Lineage2.3.4"]).
                Each value must correspond to a file named
                <strain_ids_dir>/<lineage>_strain.txt.
                When omitted all lineages in strain_ids_dir are processed.

Usage:
  snakemake --cores <N> --configfile config/config.yaml
  snakemake --cores <N> --configfile config/config.yaml --cluster "sbatch ..."
  # Run only up to SNP calling:
  snakemake --cores <N> --configfile config/config.yaml --config stop_at=step1
  # Run only up to indel calling:
  snakemake --cores <N> --configfile config/config.yaml --config stop_at=step2
  # Use a custom output directory:
  snakemake --cores <N> --configfile config/config.yaml --config outdir=results
"""

import os
import glob as _glob

configfile: "config/config.yaml"

# =============================================================================
# Global settings derived from config
# =============================================================================
OUTDIR     = config.get("outdir", "output")
STOP_AT    = config.get("stop_at", "all")
INPUT_TYPE = config.get("input_type", "auto")

# Resolve the directory where FASTQ files live.
# - fastq mode: user must explicitly provide fastq_dir (their existing files)
# - sra mode:   FASTQs are generated from SRA into <outdir>/fastq
# - auto mode:  FASTQs are looked-for / generated in <outdir>/fastq
if INPUT_TYPE == "fastq":
    if "fastq_dir" not in config:
        raise ValueError(
            "config key 'fastq_dir' must be set when input_type is 'fastq'. "
            "Set it to the directory that contains your per-sample FASTQ files."
        )
    FQ_DIR = config["fastq_dir"]
else:
    FQ_DIR = f"{OUTDIR}/fastq"

# =============================================================================
# Discover lineages and samples from strain ID files
# =============================================================================
STRAIN_IDS_DIR = config["strain_ids_dir"]

# All available lineages are derived from files named <lineage>_strain.txt
_all_lineage_names, = glob_wildcards(
    os.path.join(STRAIN_IDS_DIR, "{lineage}_strain.txt")
)
_ALL_AVAILABLE_LINEAGES = sorted(_all_lineage_names)

# If the user specified a lineage (or list of lineages), restrict to those.
_lineage_cfg = config.get("lineage", None)
if _lineage_cfg is not None:
    # Accept a single string or a list
    if isinstance(_lineage_cfg, str):
        _requested = [_lineage_cfg]
    else:
        _requested = list(_lineage_cfg)

    # Validate each requested lineage against the available strain files
    _missing = [
        lin for lin in _requested
        if not os.path.isfile(os.path.join(STRAIN_IDS_DIR, f"{lin}_strain.txt"))
    ]
    if _missing:
        raise ValueError(
            f"No strain file found for lineage(s): {_missing}. "
            f"Expected files: "
            + ", ".join(f"{STRAIN_IDS_DIR}/{m}_strain.txt" for m in _missing)
        )
    LINEAGES = sorted(_requested)
else:
    LINEAGES = _ALL_AVAILABLE_LINEAGES


def get_samples(lineage):
    """Read sample IDs from a lineage strain list file."""
    filepath = os.path.join(STRAIN_IDS_DIR, f"{lineage}_strain.txt")
    with open(filepath) as f:
        return [line.strip() for line in f if line.strip()]


ALL_SAMPLES = []
for _lin in LINEAGES:
    ALL_SAMPLES.extend(get_samples(_lin))
ALL_SAMPLES = sorted(set(ALL_SAMPLES))

# For fastq mode, filter ALL_SAMPLES to those that actually have a FASTQ file
# present in FQ_DIR, and warn about any that are missing.
if INPUT_TYPE == "fastq":
    _found = []
    _missing_fq = []
    for _s in ALL_SAMPLES:
        _pe1  = os.path.join(FQ_DIR, f"{_s}_1.fastq.gz")
        _pe2  = os.path.join(FQ_DIR, f"{_s}_2.fastq.gz")
        _se   = os.path.join(FQ_DIR, f"{_s}.fastq.gz")
        if os.path.isfile(_pe1) or os.path.isfile(_se):
            _found.append(_s)
        else:
            _missing_fq.append(_s)
    if _missing_fq:
        import sys
        print(
            f"WARNING: {len(_missing_fq)} sample(s) listed in strain files "
            f"have no matching FASTQ in {FQ_DIR!r} and will be skipped:\n"
            + "\n".join(f"  {s}" for s in _missing_fq),
            file=sys.stderr,
        )
    if not _found:
        raise ValueError(
            f"No FASTQ files found in {FQ_DIR!r} for any of the requested "
            f"samples. Check fastq_dir and the strain ID files."
        )
    ALL_SAMPLES = _found

# =============================================================================
# DR mutation selection analysis configuration (Step 7+)
# =============================================================================
DR_DIR             = config.get("dr_results_dir", "dr_results")
_DR_SNP_SOURCE     = config.get(
    "dr_convergent_snp_file",
    f"{OUTDIR}/lineage_ann/all_ann_convergent_flt.txt",
)
_DR_INDEL_FILE     = config.get("all_indel_file", "data/all_indel_100k.txt.gz")
_DR_SNP_ANNO_DIR   = config.get("snp_anno_dir", "")
_DR_INDEL_ANNO_DIR = config.get("indel_anno_dir", "")

_DR_DRUGS = [
    "RIF", "INH", "EMB", "PZA", "LFX", "MFX",
    "BDQ", "AMK", "STM", "ETO", "KAN", "CAP", "LZD",
]

# Best convergence threshold per drug (determined from training-set optimisation)
_DR_BEST_THRESHOLDS = {
    "RIF": 6, "INH": 3, "EMB": 5, "PZA": 2, "LFX": 6, "MFX": 6,
    "BDQ": 3, "AMK": 6, "STM": 5, "ETO": 4, "KAN": 6, "CAP": 6, "LZD": 5,
}

# Per-drug gene parameters used by 02/04/05 R scripts
# Keys: genes, starts, ends, lof (LoF flags for 04/05), strands (for 02/04/05)
_DR_GENE_PARAMS = {
    "RIF": {"genes": "rpoB",             "starts": "759807",                 "ends": "763325",                 "lof": "0",      "strands": "1"},
    "INH": {"genes": "inhA,katG",        "starts": "1674202,2153889",        "ends": "1675011,2156111",         "lof": "0,1",    "strands": "1,0"},
    "EMB": {"genes": "embB",             "starts": "4246514",                "ends": "4249810",                 "lof": "0",      "strands": "1"},
    "PZA": {"genes": "pncA",             "starts": "2288681",                "ends": "2289241",                 "lof": "1",      "strands": "0"},
    "LFX": {"genes": "gyrB,gyrA",        "starts": "5240,7302",              "ends": "7262,9818",               "lof": "0,0",    "strands": "1,1"},
    "MFX": {"genes": "gyrB,gyrA",        "starts": "5240,7302",              "ends": "7262,9818",               "lof": "0,0",    "strands": "1,1"},
    "BDQ": {"genes": "Rv0678,atpE,pepQ", "starts": "778990,1461045,2859300", "ends": "779487,1461290,2860418",  "lof": "1,0,1",  "strands": "1,1,0"},
    "AMK": {"genes": "rrs,eis",          "starts": "1471846,2714124",        "ends": "1473382,2715332",         "lof": "0,0",    "strands": "1,0"},
    "STM": {"genes": "gid,rpsL",         "starts": "4407528,781560",         "ends": "4408202,781934",          "lof": "1,0",    "strands": "0,1"},
    "ETO": {"genes": "inhA,ethA",        "starts": "1674202,4326004",        "ends": "1675011,4327473",         "lof": "0,1",    "strands": "1,0"},
    "KAN": {"genes": "rrs,eis",          "starts": "1471846,2714124",        "ends": "1473382,2715332",         "lof": "0,0",    "strands": "1,0"},
    "CAP": {"genes": "rrs,tlyA",         "starts": "1471846,1917940",        "ends": "1473382,1918746",         "lof": "0,1",    "strands": "1,1"},
    "LZD": {"genes": "rplC,rrl",         "starts": "800809,1473658",         "ends": "801462,1476795",          "lof": "0,0",    "strands": "1,1"},
}

# Threshold and promoter sweep parameters
_DR_THRES_SWEEP = [2, 3, 4, 5, 6]
_DR_PROM_FIXED  = 500
_DR_PROM_SWEEP  = [100, 200, 300, 400, 600, 700, 800, 900, 1000]

# All (drug, threshold, promoter) combos – threshold sweep + promoter sweep
_DR_THRES_COMBOS = [(d, t, _DR_PROM_FIXED) for d in _DR_DRUGS for t in _DR_THRES_SWEEP]
_DR_PROM_COMBOS  = [(d, _DR_BEST_THRESHOLDS[d], p) for d in _DR_DRUGS for p in _DR_PROM_SWEEP]
_DR_ALL_COMBOS   = _DR_THRES_COMBOS + _DR_PROM_COMBOS

# Output targets for each DR analysis stage
_DR_PREP_TARGETS = (
    expand(f"{DR_DIR}/{{drug}}/id/train_70.txt",                      drug=_DR_DRUGS) +
    expand(f"{DR_DIR}/{{drug}}/denovo_snp_2.txt",                     drug=_DR_DRUGS) +
    expand(f"{DR_DIR}/{{drug}}/denovo_EvoResist_initial_list.txt",    drug=_DR_DRUGS)
)

_DR_SWEEP_TARGETS = _DR_PREP_TARGETS + [
    f"{DR_DIR}/{d}/Threshold_{t}_Promoter_{p}_list2.tsv"
    for d, t, p in _DR_ALL_COMBOS
] + [
    f"{DR_DIR}/{d}/Threshold_{t}_Promoter_{p}/train_list2/overall_metrics.tsv"
    for d, t, p in _DR_ALL_COMBOS
] + [
    f"{DR_DIR}/{d}/Threshold_{t}_Promoter_{p}/test_list2/overall_metrics.tsv"
    for d, t, p in _DR_ALL_COMBOS
]

# =============================================================================
# Map each stop_at value to its corresponding output files
# =============================================================================
_SIMULATION_OUTPUTS = [
    f"{OUTDIR}/simulation/simulated_mutations_raw_GTR_Gamma.csv",
    f"{OUTDIR}/simulation/null_mutation_df_GTR.csv",
    f"{OUTDIR}/simulation/expected_null_distribution_GTR_Gamma.csv",
]

_INDEL_OUTPUTS = expand(f"{OUTDIR}/indel/{{sample}}.INDEL.pass.vcf.gz", sample=ALL_SAMPLES)

_STOP_STEP_MAP = {
    "step1":              expand(f"{OUTDIR}/snv/{{sample}}.snp", sample=ALL_SAMPLES),
    "snp_calling":        expand(f"{OUTDIR}/snv/{{sample}}.snp", sample=ALL_SAMPLES),
    "step2":              _INDEL_OUTPUTS,
    "indel_calling":      _INDEL_OUTPUTS,
    "step3":              expand(f"{OUTDIR}/lineage_tree/{{lineage}}_btp.treefile", lineage=LINEAGES),
    "build_tree":         expand(f"{OUTDIR}/lineage_tree/{{lineage}}_btp.treefile", lineage=LINEAGES),
    "step4":              expand(f"{OUTDIR}/lineage_ann/{{lineage}}.ann", lineage=LINEAGES),
    "branch_mutations":   expand(f"{OUTDIR}/lineage_ann/{{lineage}}.ann", lineage=LINEAGES),
    "step5":              [f"{OUTDIR}/lineage_ann/ancestor.ann"],
    "ancestor_mutations": [f"{OUTDIR}/lineage_ann/ancestor.ann"],
    "step6":              [f"{OUTDIR}/lineage_ann/all_ann_convergent_flt.txt"],
    "merge_annotations":  [f"{OUTDIR}/lineage_ann/all_ann_convergent_flt.txt"],
    "step7":              _SIMULATION_OUTPUTS,
    "simulation":         _SIMULATION_OUTPUTS,
    "all":                _SIMULATION_OUTPUTS,
    # DR mutation selection analysis (Step 8+)
    "step8":              _DR_PREP_TARGETS,
    "dr_prep":            _DR_PREP_TARGETS,
    "step9":              _DR_SWEEP_TARGETS,
    "dr_selection":       _DR_SWEEP_TARGETS,
}

if STOP_AT not in _STOP_STEP_MAP:
    raise ValueError(
        f"Invalid stop_at value '{STOP_AT}'. "
        f"Valid options: {sorted(set(_STOP_STEP_MAP.keys()))}"
    )

# Validate that DR analysis prerequisites are set when a DR stop is requested
if STOP_AT in ("step8", "dr_prep", "step9", "dr_selection"):
    if not _DR_SNP_ANNO_DIR:
        raise ValueError(
            "config key 'snp_anno_dir' must be set when stop_at is a DR analysis step. "
            "Point it to the directory containing per-sample SNP annotation files "
            "named {sample}.ano (columns: position, ref, alt)."
        )
    if not _DR_INDEL_ANNO_DIR:
        raise ValueError(
            "config key 'indel_anno_dir' must be set when stop_at is a DR analysis step. "
            "Point it to the directory containing per-sample indel annotation files "
            "named {sample}.indel.ano (columns: position, ref, alt)."
        )

FINAL_TARGETS = _STOP_STEP_MAP[STOP_AT]

# =============================================================================
# Final target
# =============================================================================
rule all:
    input: FINAL_TARGETS


# =============================================================================
# Step 1: SNP Calling (per sample)
# =============================================================================
# Processes FASTQ files through quality trimming (sickle), alignment (bwa),
# variant calling (VarScan), and filtering to produce SNP, CFA, and forup files.

rule snp_calling:
    input:
        ref=config["reference"],
        ref_fai=config["reference"] + ".fai",
        low_ebr=config["low_ebr_file"],
    output:
        snp=f"{OUTDIR}/snv/{{sample}}.snp",
        cfa=f"{OUTDIR}/cfa/{{sample}}.cfa",
        forup=f"{OUTDIR}/forup/{{sample}}.forup",
        bam=f"{OUTDIR}/bam/{{sample}}.sort.bam",
    params:
        sample="{sample}",
        fq_dir=FQ_DIR,
        bam_dir=f"{OUTDIR}/bam",
        snv_dir=f"{OUTDIR}/snv",
        cfa_dir=f"{OUTDIR}/cfa",
        forup_dir=f"{OUTDIR}/forup",
        varscan=config["varscan_jar"],
        sra_dir=config["sra_dir"],
        input_type=INPUT_TYPE,
    threads: config["threads_per_sample"]
    resources:
        mem_mb=12000,
    shell:
        r"""
        set -euo pipefail
        sample={params.sample}
        ref={input.ref}
        fq_dir={params.fq_dir}
        bam_dir={params.bam_dir}
        snv={params.snv_dir}
        cfa_dir={params.cfa_dir}
        forup_dir={params.forup_dir}
        varscan={params.varscan}
        sra_dir={params.sra_dir}
        low_ebr={input.low_ebr}
        nthreads={threads}
        input_type={params.input_type}

        mkdir -p "$bam_dir" "$snv" "$cfa_dir" "$forup_dir"
        # Only create fq_dir when it's an output location (sra / auto modes)
        if [ "$input_type" != "fastq" ]; then
            mkdir -p "$fq_dir"
        fi

        # -- Locate FASTQ files --
        fq1="${{fq_dir}}/${{sample}}_1.fastq.gz"
        fq2="${{fq_dir}}/${{sample}}_2.fastq.gz"
        sfq="${{fq_dir}}/${{sample}}.fastq.gz"

        # -- Resolve input based on input_type --
        if [ "$input_type" = "sra" ]; then
            # SRA mode: convert SRA to FASTQ unconditionally
            if [ -f "${{sra_dir}}/${{sample}}/${{sample}}.sra" ]; then
                fastq-dump --split-3 --gzip -O "$fq_dir" \
                    "${{sra_dir}}/${{sample}}/${{sample}}.sra"
            elif [ -f "${{sra_dir}}/${{sample}}/${{sample}}.sralite" ]; then
                fastq-dump --split-3 --gzip -O "$fq_dir" \
                    "${{sra_dir}}/${{sample}}/${{sample}}.sralite"
            else
                echo "Error: No SRA file found for $sample in $sra_dir" >&2
                exit 1
            fi
        elif [ "$input_type" = "fastq" ]; then
            # FASTQ mode: expect files to already exist
            if [ ! -f "$fq1" ] && [ ! -f "$sfq" ]; then
                echo "Error: No FASTQ file found for $sample. Expected $fq1 or $sfq" >&2
                exit 1
            fi
        else
            # Auto mode: use FASTQ if present, otherwise fall back to SRA
            if [ ! -f "$fq1" ] && [ ! -f "$sfq" ]; then
                if [ -f "${{sra_dir}}/${{sample}}/${{sample}}.sra" ]; then
                    fastq-dump --split-3 --gzip -O "$fq_dir" \
                        "${{sra_dir}}/${{sample}}/${{sample}}.sra"
                elif [ -f "${{sra_dir}}/${{sample}}/${{sample}}.sralite" ]; then
                    fastq-dump --split-3 --gzip -O "$fq_dir" \
                        "${{sra_dir}}/${{sample}}/${{sample}}.sralite"
                else
                    echo "Error: No FASTQ or SRA file found for $sample" >&2
                    exit 1
                fi
            fi
        fi

        # -- Define intermediate file paths --
        sortbam="${{bam_dir}}/${{sample}}.sort.bam"
        pileup="${{snv}}/${{sample}}.pileup"
        var="${{snv}}/${{sample}}.varscan"
        cns="${{snv}}/${{sample}}.cns"
        ppe="${{snv}}/${{sample}}.ppe"
        format="${{snv}}/${{sample}}.for"
        forup="${{forup_dir}}/${{sample}}.forup"
        fix="${{snv}}/${{sample}}.fix"
        snp="${{snv}}/${{sample}}.snp"
        cfa="${{cfa_dir}}/${{sample}}.cfa"

        if [ -f "$fq2" ]; then
            # === Paired-end processing ===
            fq1_tr="${{fq_dir}}/${{sample}}_tr_1.fq"
            fq2_tr="${{fq_dir}}/${{sample}}_tr_2.fq"
            fq3_tr="${{fq_dir}}/${{sample}}_tr_S.fq"
            sai1="${{fq_dir}}/${{sample}}_R1.sai"
            sai2="${{fq_dir}}/${{sample}}_R2.sai"
            sai3="${{fq_dir}}/${{sample}}_S.sai"
            samp="${{bam_dir}}/${{sample}}.paired.sam"
            sams="${{bam_dir}}/${{sample}}.single.sam"
            bamp="${{bam_dir}}/${{sample}}.paired.bam"
            bams="${{bam_dir}}/${{sample}}.single.bam"
            bamm="${{bam_dir}}/${{sample}}.merge.bam"

            sickle pe -l 35 -f "$fq1" -r "$fq2" -t sanger \
                -o "$fq1_tr" -p "$fq2_tr" -s "$fq3_tr"
            bwa aln -R 1 -t "$nthreads" "$ref" "$fq1_tr" > "$sai1"
            bwa aln -R 1 -t "$nthreads" "$ref" "$fq2_tr" > "$sai2"
            bwa aln -R 1 "$ref" "$fq3_tr" > "$sai3"
            bwa sampe -a 1000 -n 1 -N 1 "$ref" \
                "$sai1" "$sai2" "$fq1_tr" "$fq2_tr" > "$samp"
            bwa samse -n 1 "$ref" "$sai3" "$fq3_tr" > "$sams"
            samtools view -bhSt "${{ref}}.fai" "$samp" -o "$bamp"
            samtools view -bhSt "${{ref}}.fai" "$sams" -o "$bams"
            samtools merge "$bamm" "$bamp" "$bams"
            samtools sort "$bamm" -o "$sortbam"
            samtools index "$sortbam"
            samtools mpileup -q 40 -Q 30 -ABOf "$ref" "$sortbam" > "$pileup"

            java -jar "$varscan" mpileup2snp "$pileup" \
                --min-coverage 10 --min-reads2 4 --min-avg-qual 30 \
                --min-var-freq 0.01 --min-freq-for-hom 0.9 \
                --p-value 99e-02 --strand-filter 0 > "$var"
            java -jar "$varscan" mpileup2cns "$pileup" \
                --min-coverage 10 --min-avg-qual 30 --min-var-freq 0.9 \
                --min-reads2 4 --strand-filter 0 > "$cns"

            python scripts/remove_low_ebr.py "$low_ebr" "$var" > "$ppe"
            perl scripts/varscan_work_flow/1_format_trans.pl "$ppe" > "$format"
            perl scripts/varscan_work_flow/2_fix_extract.pl "$format" > "$fix"
            perl scripts/varscan_work_flow/3.1_mix_pileup_merge.pl \
                "$format" "$pileup" > "$forup"
            cut -f2-4 "$fix" > "$snp"
            perl scripts/single_colony/1st_loci_recall_cns.pl "$cns" > "$cfa"

            rm -f "$fq1_tr" "$fq2_tr" "$fq3_tr" \
                  "$sai1" "$sai2" "$sai3" "$samp" "$sams" \
                  "$bamp" "$bams" "$bamm" "$pileup" "$ppe" "$format" "$fix"
            rm -f "${{snv}}/${{sample}}.cns.gz"
            gzip "$cns"
        else
            # === Single-end processing ===
            fq_tr="${{fq_dir}}/${{sample}}_tr.fq"
            sai="${{fq_dir}}/${{sample}}.sai"
            samf="${{bam_dir}}/${{sample}}.sam"
            bamf="${{bam_dir}}/${{sample}}.bam"

            sickle se -f "$sfq" -t sanger -o "$fq_tr"
            bwa aln -t "$nthreads" -R 1 "$ref" "$fq_tr" > "$sai"
            bwa samse "$ref" "$sai" "$fq_tr" > "$samf"
            samtools view -bhSt "${{ref}}.fai" "$samf" -o "$bamf"
            samtools sort "$bamf" -o "$sortbam"
            samtools index "$sortbam"
            samtools mpileup -q 40 -Q 30 -ABOf "$ref" "$sortbam" > "$pileup"

            java -jar "$varscan" mpileup2snp "$pileup" \
                --min-coverage 10 --min-reads2 4 --min-avg-qual 30 \
                --min-var-freq 0.01 --min-freq-for-hom 0.9 \
                --p-value 99e-02 --strand-filter 0 > "$var"
            java -jar "$varscan" mpileup2cns "$pileup" \
                --min-coverage 10 --min-avg-qual 30 --min-var-freq 0.9 \
                --min-reads2 4 --strand-filter 0 > "$cns"

            python scripts/remove_low_ebr.py "$low_ebr" "$var" > "$ppe"
            perl scripts/varscan_work_flow/1_format_trans.pl "$ppe" > "$format"
            perl scripts/varscan_work_flow/2_fix_extract.pl "$format" > "$fix"
            perl scripts/varscan_work_flow/3.1_mix_pileup_merge.pl \
                "$format" "$pileup" > "$forup"
            cut -f2-4 "$fix" > "$snp"
            perl scripts/single_colony/1st_loci_recall_cns.pl "$cns" > "$cfa"

            rm -f "$fq_tr" "$sai" "$samf" "$bamf" \
                  "$pileup" "$ppe" "$format" "$fix"
            rm -f "${{snv}}/${{sample}}.cns.gz"
            gzip "$cns"
        fi
        echo "SNP calling complete for $sample."
        """


# =============================================================================
# Step 2: INDEL Calling (per sample)
# =============================================================================
# Calls INDELs from sorted BAM files using GATK3 local realignment and
# UnifiedGenotyper, followed by bcftools normalization and QC filtering.
# Produces per-sample INDEL VCF files (pass-filtered) in {OUTDIR}/indel/.

rule indel_calling:
    input:
        bam=f"{OUTDIR}/bam/{{sample}}.sort.bam",
        ref=config["reference"],
    output:
        vcf_pass=f"{OUTDIR}/indel/{{sample}}.INDEL.pass.vcf.gz",
        vcf_pass_tbi=f"{OUTDIR}/indel/{{sample}}.INDEL.pass.vcf.gz.tbi",
    params:
        sample="{sample}",
        indel_dir=f"{OUTDIR}/indel",
        picard=config["picard_jar"],
        gatk=config["gatk_jar"],
    resources:
        mem_mb=8000,
    shell:
        r"""
        set -euo pipefail
        sample={params.sample}
        indel_dir={params.indel_dir}
        ref={input.ref}
        bam={input.bam}
        picard={params.picard}
        gatk={params.gatk}

        mkdir -p "$indel_dir"

        # QC thresholds
        MAX_INDEL_LEN=50
        MIN_DP=10
        MIN_AF=0.75

        # Intermediate BAM files
        bam_filt="${{indel_dir}}/${{sample}}.filt.bam"
        bam_rg="${{indel_dir}}/${{sample}}.RG.bam"
        intervals="${{indel_dir}}/${{sample}}.intervals"
        bam_rl="${{indel_dir}}/${{sample}}.RL.bam"

        # VCF output files
        vcf_raw="${{indel_dir}}/${{sample}}.INDEL.raw.vcf"
        vcf_norm="${{indel_dir}}/${{sample}}.INDEL.norm.vcf.gz"
        vcf_tags="${{indel_dir}}/${{sample}}.INDEL.norm.tags.vcf.gz"
        vcf_pass="{output.vcf_pass}"

        ###########################################################################
        # 1) Filter BAM: remove unmapped reads (-F 4) and MAPQ=0 reads (-q 1)
        ###########################################################################
        samtools view -bF 4 -q 1 "$bam" > "$bam_filt"
        samtools index "$bam_filt"

        ###########################################################################
        # 2) Add or replace read groups (required by GATK tools)
        ###########################################################################
        java -jar "$picard" AddOrReplaceReadGroups \
            I="$bam_filt" O="$bam_rg" \
            RGID=4 RGLB=lib1 RGPL=illumina RGPU=unit1 RGSM=20
        samtools index "$bam_rg"

        ###########################################################################
        # 3) GATK3: identify realignment targets around indels
        ###########################################################################
        java -Xmx2g -jar "$gatk" \
            -I "$bam_rg" -R "$ref" \
            -T RealignerTargetCreator -maxInterval 20000 \
            -o "$intervals"

        ###########################################################################
        # 4) GATK3: local realignment around indels
        ###########################################################################
        java -Xmx4g -jar "$gatk" \
            -I "$bam_rg" -R "$ref" \
            -T IndelRealigner -targetIntervals "$intervals" \
            -o "$bam_rl"
        samtools index "$bam_rl"

        ###########################################################################
        # 5) GATK3: call INDELs (UnifiedGenotyper)
        ###########################################################################
        java -Xmx4g -jar "$gatk" \
            -T UnifiedGenotyper -nt 1 \
            -R "$ref" -I "$bam_rl" \
            -glm INDEL -o "$vcf_raw"

        ###########################################################################
        # 6) Normalize INDEL representation (left-align; split multiallelic)
        ###########################################################################
        bcftools norm -f "$ref" -m -both "$vcf_raw" -Oz -o "$vcf_norm"
        bcftools index -t "$vcf_norm"

        ###########################################################################
        # 7) Add tags AF/AC/AN (computed from genotype fields)
        ###########################################################################
        bcftools +fill-tags "$vcf_norm" -Oz -o "$vcf_tags" -- -t AF,AC,AN
        bcftools index -t "$vcf_tags"

        ###########################################################################
        # 8) QC filter:
        #    - keep only INDELs
        #    - indel length <= 50 bp
        #    - FORMAT/DP >= 10  (per-sample depth)
        #    - INFO/AF >= 0.75  (fixed/high-confidence variants)
        ###########################################################################
        bcftools view -i "TYPE='indel' && abs(strlen(REF)-strlen(ALT))<=${{MAX_INDEL_LEN}}" \
            "$vcf_tags" -Ou \
            | bcftools filter -i "FORMAT/DP>=${{MIN_DP}} && INFO/AF>=${{MIN_AF}}" \
            -Oz -o "$vcf_pass"
        bcftools index -t "$vcf_pass"

        ###########################################################################
        # 9) Cleanup intermediate BAM files
        ###########################################################################
        rm -f "$bam_filt" "$bam_filt.bai" \
              "$bam_rg" "$bam_rg.bai" \
              "$intervals" \
              "$bam_rl" "$bam_rl.bai"

        echo "INDEL calling complete for $sample."
        """


# =============================================================================
# Step 3: Build Phylogenetic Tree (per lineage/sublineage)
# =============================================================================
# Collects SNP positions across all samples in a lineage, generates a
# concatenated alignment, and builds a maximum-likelihood tree with IQ-TREE.

rule build_tree:
    input:
        snps=lambda wc: expand(f"{OUTDIR}/snv/{{sample}}.snp",
                               sample=get_samples(wc.lineage)),
        cfas=lambda wc: expand(f"{OUTDIR}/cfa/{{sample}}.cfa",
                               sample=get_samples(wc.lineage)),
        strain_list=os.path.join(config["strain_ids_dir"], "{lineage}_strain.txt"),
        anc_cfa=config["ancestor_concat_fasta"],
        ref=config["reference"],
    output:
        tree=f"{OUTDIR}/lineage_tree/{{lineage}}_btp.treefile",
        state=f"{OUTDIR}/lineage_tree/{{lineage}}_btp.state",
    params:
        lineage="{lineage}",
        tree_dir=f"{OUTDIR}/lineage_tree",
        cfa_dir=f"{OUTDIR}/cfa",
        snv_dir=f"{OUTDIR}/snv",
    threads: config["threads_tree"]
    resources:
        mem_mb=200000,
    shell:
        r"""
        set -euo pipefail
        lineage={params.lineage}
        tree_dir={params.tree_dir}
        cfa_dir={params.cfa_dir}
        snv_dir={params.snv_dir}
        mkdir -p "$tree_dir"

        # Generate position file from all samples' SNPs
        pos_file="${{tree_dir}}/${{lineage}}.pos"
        snp_strain_file="${{tree_dir}}/${{lineage}}_snp_strain.txt"
        cfa_strain_file="${{tree_dir}}/${{lineage}}_cfa_strain.txt"
        > "$pos_file"
        > "$snp_strain_file"
        > "$cfa_strain_file"

        while IFS= read -r strain; do
            snp_file="${{snv_dir}}/${{strain}}.snp"
            if [ -f "$snp_file" ]; then
                awk '{{print $1}}' "$snp_file" >> "$pos_file"
                echo "$strain" >> "$snp_strain_file"
            else
                echo "Warning: ${{strain}} snp file not found, skipping." >&2
            fi
        done < {input.strain_list}

        sort -nu "$pos_file" -o "$pos_file"

        # Generate CFA strain list
        while IFS= read -r strain; do
            cfa_file="${{cfa_dir}}/${{strain}}.cfa"
            if [ -f "$cfa_file" ]; then
                echo "$cfa_file" >> "$cfa_strain_file"
            else
                echo "Warning: ${{strain}} cfa file not found, skipping." >&2
            fi
        done < "$snp_strain_file"
        echo "{input.anc_cfa}" >> "$cfa_strain_file"

        # Generate concatenated alignment
        python scripts/snp2cfa.py "$pos_file" "$cfa_strain_file"

        # Build tree with IQ-TREE
        cfa="${{tree_dir}}/${{lineage}}.fa"
        if [ ! -f "$cfa" ]; then
            echo "Error: Expected CFA file $cfa not found!" >&2
            exit 1
        fi

        iqtree2 -s "$cfa" \
            -m GTR+F+R4 \
            --seqtype DNA \
            --prefix "${{tree_dir}}/${{lineage}}_btp" \
            --mem {resources.mem_mb}M \
            -T {threads} \
            --ancestral \
            -o tb \
            -af fasta \
            -bb 1000 \
            -alrt 1000
        """


# =============================================================================
# Step 4: Branch Mutation Extraction (within lineage/sublineage)
# =============================================================================
# Extracts mutations per branch/node from the phylogenetic tree and annotates
# them using the MTBC translation scripts.

rule branch_mutations:
    input:
        tree=f"{OUTDIR}/lineage_tree/{{lineage}}_btp.treefile",
        state=f"{OUTDIR}/lineage_tree/{{lineage}}_btp.state",
        low_ebr=config["low_ebr_file"],
        ref=config["reference"],
    output:
        ann=f"{OUTDIR}/lineage_ann/{{lineage}}.ann",
    params:
        lineage="{lineage}",
        tree_dir=f"{OUTDIR}/lineage_tree",
        ann_dir=f"{OUTDIR}/lineage_ann",
        cfa_dir=f"{OUTDIR}/lineage_cfa",
    threads: config["threads_annotation"]
    resources:
        mem_mb=20000,
    shell:
        r"""
        set -euo pipefail
        lineage={params.lineage}
        tree_dir={params.tree_dir}
        ann_dir={params.ann_dir}
        cfa_dir={params.cfa_dir}
        low_ebr={input.low_ebr}
        mkdir -p "$ann_dir"

        cd "$tree_dir"

        # Extract ancestral sequences for each node
        perl ../../scripts/nodes_base_locus_iqtree.pl \
            ${{lineage}}_btp.treefile \
            ../../${{cfa_dir}}/${{lineage}}_delete.pos \
            ${{lineage}}_btp.state \
            ${{lineage}}.fa \
            ${{lineage}}.db \
            ${{lineage}}_homoplasy.txt

        # Extract mutations from the database file
        perl ../../scripts/db2mutation.pl ${{lineage}}.db \
            > ${{lineage}}_db_mutation.txt

        # Annotate with leaf counts per node
        python ../../scripts/node_leafs.py \
            ${{lineage}}_btp.treefile \
            ${{lineage}}_db_mutation.txt \
            > ${{lineage}}_db_mutation2.txt

        # Remove ancestral (lineage-defining) mutations
        python ../../scripts/filter_lineage_defining.py \
            ${{lineage}}_db_mutation2.txt \
            ${{lineage}} \
            > ${{lineage}}_db_mutation2_rmanc.txt

        # Generate per-node SNP files
        python ../../scripts/getrefbase_per_node.py \
            ${{lineage}}_db_mutation2_rmanc.txt ${{lineage}} \
            ../../{input.ref}

        # Annotate each node's mutations in parallel
        find "${{lineage}}" -maxdepth 1 -name "*snp" -print0 | \
            parallel -0 -j {threads} '
                python ../../scripts/remove_low_ebr.py \
                    "../../'"$low_ebr"'" "{{}}" \
                    > "{{= s/\.snp$// =}}_rle.snp" &&
                perl ../../scripts/mtbc_translate/0_MTBC_Annotation_mtbc_4411532_corrected.pl \
                    "{{= s/\.snp$// =}}_rle.snp" \
                    > "{{= s/\.snp$// =}}.ann" &&
                sed -i "/^$/d" "{{= s/\.snp$// =}}.ann"
            '

        cat "${{lineage}}"/*.ann > "../../${{ann_dir}}/${{lineage}}.ann"

        # Clean up temporary files
        rm -f ${{lineage}}_db_mutation.txt ${{lineage}}_db_mutation2.txt
        rm -rf "${{lineage}}"

        cd ../..
        """


# =============================================================================
# Step 5: Ancestor Mutation Extraction (prior to lineage/sublineage diversification)
# =============================================================================
# Extracts and annotates mutations from ancestor nodes shared across lineages.

rule ancestor_mutations:
    input:
        mutation_file=config["ancestor_mutation_file"],
        low_ebr=config["low_ebr_file"],
        ref=config["reference"],
    output:
        ann=f"{OUTDIR}/lineage_ann/ancestor.ann",
    params:
        ann_dir=f"{OUTDIR}/lineage_ann",
        tmp_dir=f"{OUTDIR}/ancestor_tmp",
    threads: config["threads_annotation"]
    resources:
        mem_mb=20000,
    shell:
        r"""
        set -euo pipefail
        ann_dir={params.ann_dir}
        tmp_dir={params.tmp_dir}
        low_ebr={input.low_ebr}
        mkdir -p "$ann_dir"

        # Generate per-node SNP files for ancestor nodes
        python scripts/getrefbase_per_node.py \
            {input.mutation_file} "$tmp_dir" {input.ref}

        # Annotate each node's mutations in parallel
        find "$tmp_dir" -maxdepth 1 -name "*snp" -print0 | \
            parallel -0 -j {threads} '
                python scripts/remove_low_ebr.py "'"$low_ebr"'" "{{}}" \
                    > "{{= s/\.snp$// =}}_rle.snp" &&
                perl scripts/mtbc_translate/0_MTBC_Annotation_mtbc_4411532_corrected.pl \
                    "{{= s/\.snp$// =}}_rle.snp" \
                    > "{{= s/\.snp$// =}}.ann" &&
                sed -i "/^$/d" "{{= s/\.snp$// =}}.ann"
            '

        cat "$tmp_dir"/*.ann > "$ann_dir/ancestor.ann"

        # Clean up
        rm -rf "$tmp_dir"
        """


# =============================================================================
# Step 6: Merge Annotations and Count Convergent Mutations by Codon
# =============================================================================
# Merges all lineage and ancestor annotations, then counts convergent events.

rule merge_annotations:
    input:
        lineage_anns=expand(f"{OUTDIR}/lineage_ann/{{lineage}}.ann",
                            lineage=LINEAGES),
        ancestor_ann=f"{OUTDIR}/lineage_ann/ancestor.ann",
    output:
        all_ann=f"{OUTDIR}/lineage_ann/all_ann.txt",
    shell:
        r"""
        cat {input.lineage_anns} > {output.all_ann}
        cat {input.ancestor_ann} >> {output.all_ann}
        """


rule stat_convergent:
    input:
        all_ann=f"{OUTDIR}/lineage_ann/all_ann.txt",
    output:
        convergent=f"{OUTDIR}/lineage_ann/all_ann_convergent.txt",
    params:
        ann_dir=f"{OUTDIR}/lineage_ann",
    shell:
        r"""
        ann_dir={params.ann_dir}
        cd "$ann_dir"
        Rscript ../../scripts/stat_convergent.R
        """


rule filter_convergent:
    input:
        convergent=f"{OUTDIR}/lineage_ann/all_ann_convergent.txt",
        snp_freq=config["snp_freq_file"],
        repeat_region=config["repeat_region_file"],
        mobile_element=config["mobile_element_file"],
    output:
        filtered=f"{OUTDIR}/lineage_ann/all_ann_convergent_flt.txt",
    params:
        ann_dir=f"{OUTDIR}/lineage_ann",
    shell:
        r"""
        ann_dir={params.ann_dir}
        cd "$ann_dir"
        Rscript ../../scripts/filter_low_freq_pos.R \
            all_ann_convergent.txt all_ann_convergent_flt.txt
        """


# =============================================================================
# Step 7: GTR+Gamma Simulation of Mutations Under a Null Distribution
# =============================================================================
# Simulates convergent mutation counts under a GTR+Gamma model to generate
# the null distribution for statistical testing.

rule simulation:
    input:
        filtered=f"{OUTDIR}/lineage_ann/all_ann_convergent_flt.txt",
        ref=config["reference"],
    output:
        raw_sim=f"{OUTDIR}/simulation/simulated_mutations_raw_GTR_Gamma.csv",
        null_dist=f"{OUTDIR}/simulation/null_mutation_df_GTR.csv",
        expected=f"{OUTDIR}/simulation/expected_null_distribution_GTR_Gamma.csv",
    params:
        sim_dir=f"{OUTDIR}/simulation",
    shell:
        r"""
        mkdir -p {params.sim_dir}
        cd {params.sim_dir}
        #optional, if you want to get GTR probs using your own data, you can run the following command 
        #and replace the output GTR probs in the simulation_GTR_gamma.py script
        #python ../../scripts/fourfold_dgr_rate.py ../../{OUTDIR}/lineage_ann/
        
        python ../../scripts/simulation_GTR_gamma.py
        """


# =============================================================================
# Step 8: DR – Stratify 70/30 Train-Test Split per Drug
# =============================================================================
# Stratified 70/30 train-test split of the per-drug sample lists for
# sensitivity analysis to identify drug-specific convergent thresholds.

rule dr_train_test_split:
    input:
        sample_list="data/{drug}_sample_list.txt",
    output:
        train=f"{DR_DIR}/{{drug}}/id/train_70.txt",
        test=f"{DR_DIR}/{{drug}}/id/test_30.txt",
    params:
        outdir=lambda wc: f"{DR_DIR}/{wc.drug}/id",
    shell:
        r"""
        set -euo pipefail
        Rscript scripts/dr_mutation_selection/01-train_test_split.R \
            {input.sample_list} {params.outdir}
        """


# =============================================================================
# Step 9: DR – Filter Convergent Variants to Drug Resistance Genes and Promoter Regions
# =============================================================================
# Applies drug-specific genomic region filters to include convergent SNPs or
# indels from drug resistance genes and promoter regions.

rule dr_filter_variants:
    input:
        snp_file=_DR_SNP_SOURCE,
        indel_file=_DR_INDEL_FILE,
    output:
        snp=f"{DR_DIR}/{{drug}}/denovo_snp_2.txt",
        indel=f"{DR_DIR}/{{drug}}/denovo_indel_2.txt",
    params:
        drug="{drug}",
        outdir=lambda wc: f"{DR_DIR}/{wc.drug}",
    shell:
        r"""
        set -euo pipefail
        mkdir -p {params.outdir}

        # Filter SNPs: all positions in each (possibly range-encoded) entry must
        # fall within the drug-relevant genomic region.
        # NOTE: SNP windows intentionally extend ~1000 bp upstream to capture
        # promoter-region variants; indel filters below use exact gene coordinates.
        # Boundaries match the filter_snp_allpos_inside function in 00-pipeline.sh.
        snp_file={input.snp_file}
        [[ "$snp_file" == *.gz ]] && snp_cat="zcat" || snp_cat="cat"
        $snp_cat "$snp_file" | awk -F' ' -v DRUG='{params.drug}' '
        $4 ~ /^[0-9]+(-[0-9]+)*$/ {{
            n = split($4, a, "-")
            keep = 1
            for (i = 1; i <= n; i++) {{
                pos = a[i]+0
                ok = 0
                if      (DRUG == "RIF") {{ ok = (pos >= 758807  && pos <= 763325) }}
                else if (DRUG == "INH") {{ ok = ((pos >= 1672440 && pos <= 1675011) || (pos >= 2153889 && pos <= 2157111)) }}
                else if (DRUG == "EMB") {{ ok = (pos >= 4245514  && pos <= 4249810) }}
                else if (DRUG == "PZA") {{ ok = (pos >= 2288681  && pos <= 2290241) }}
                else if (DRUG == "LFX" || DRUG == "MFX") {{ ok = (pos >= 4240 && pos <= 9818) }}
                else if (DRUG == "BDQ") {{ ok = ((pos >= 777990  && pos <= 779487)  || (pos >= 1460045 && pos <= 1461290) || (pos >= 2858300 && pos <= 2861418)) }}
                else if (DRUG == "AMK") {{ ok = ((pos >= 1470846 && pos <= 1473382) || (pos >= 2715332 && pos <= 2716332)) }}
                else if (DRUG == "STM") {{ ok = ((pos >= 780560  && pos <= 781934)  || (pos >= 4407528 && pos <= 4409202) || pos == 1472359 || pos == 1472362) }}
                else if (DRUG == "ETO") {{ ok = ((pos >= 4326004 && pos <= 4328473) || (pos >= 1672440 && pos <= 1675011)) }}
                else if (DRUG == "KAN") {{ ok = ((pos >= 1470846 && pos <= 1473382) || (pos >= 2714124 && pos <= 2716332)) }}
                else if (DRUG == "CAP") {{ ok = ((pos >= 1470846 && pos <= 1473382) || (pos >= 1916940 && pos <= 1918746)) }}
                else if (DRUG == "LZD") {{ ok = ((pos >= 799809  && pos <= 801462)  || (pos >= 1472658 && pos <= 1476795)) }}
                if (!ok) {{ keep = 0; break }}
            }}
            if (keep) print
        }}' > {output.snp}

        # Filter indels by position (column 2, tab-separated)
        indel_file={input.indel_file}
        [[ "$indel_file" == *.gz ]] && indel_cat="zcat" || indel_cat="cat"
        $indel_cat "$indel_file" | awk -F'\t' -v DRUG='{params.drug}' '{{
            pass = 0
            if      (DRUG == "RIF") {{ pass = ($2 >= 759807  && $2 <= 763325) }}
            else if (DRUG == "INH") {{ pass = (($2 >= 1674202 && $2 <= 1675011) || ($2 >= 2153889 && $2 <= 2156111)) }}
            else if (DRUG == "EMB") {{ pass = ($2 >= 4246514  && $2 <= 4249810) }}
            else if (DRUG == "PZA") {{ pass = ($2 >= 2288681  && $2 <= 2289241) }}
            else if (DRUG == "LFX" || DRUG == "MFX") {{ pass = (($2 >= 5240 && $2 <= 7262) || ($2 >= 7302 && $2 <= 9818)) }}
            else if (DRUG == "BDQ") {{ pass = (($2 >= 778990  && $2 <= 779487)  || ($2 >= 1461045 && $2 <= 1461290) || ($2 >= 2859300 && $2 <= 2860418)) }}
            else if (DRUG == "AMK") {{ pass = (($2 >= 1471846 && $2 <= 1473382) || ($2 >= 2715332 && $2 <= 2716332)) }}
            else if (DRUG == "STM") {{ pass = (($2 >= 781560  && $2 <= 781934)  || ($2 >= 4407528 && $2 <= 4408202)) }}
            else if (DRUG == "ETO") {{ pass = (($2 >= 4326004 && $2 <= 4327473) || ($2 >= 1674202 && $2 <= 1675011)) }}
            else if (DRUG == "KAN") {{ pass = (($2 >= 1471846 && $2 <= 1473382) || ($2 >= 2714124 && $2 <= 2716332)) }}
            else if (DRUG == "CAP") {{ pass = (($2 >= 1471846 && $2 <= 1473382) || ($2 >= 1917940 && $2 <= 1918746)) }}
            else if (DRUG == "LZD") {{ pass = (($2 >= 800809  && $2 <= 801462)  || ($2 >= 1473658 && $2 <= 1476795)) }}
            if (pass) print
        }}' > {output.indel}
        """


# =============================================================================
# Step 10: DR – Build Annotated Initial Candidate Variant List (per drug)
# =============================================================================
# Merges convergent SNPs and indels with WHO catalogue annotations, assigns
# gene-body / promoter labels, and aggregates convergence event counts.
# Requires per-drug WHO catalogue files placed at:
#   <dr_results_dir>/<drug>/WHO_list/WHO_list_allGroup.txt

rule dr_initial_list:
    input:
        snp=f"{DR_DIR}/{{drug}}/denovo_snp_2.txt",
        indel=f"{DR_DIR}/{{drug}}/denovo_indel_2.txt",
        who=f"{DR_DIR}/{{drug}}/WHO_list/WHO_list_allGroup.txt",
    output:
        f"{DR_DIR}/{{drug}}/denovo_EvoResist_initial_list.txt",
    params:
        drug="{drug}",
        genes=lambda wc:   _DR_GENE_PARAMS[wc.drug]["genes"],
        starts=lambda wc:  _DR_GENE_PARAMS[wc.drug]["starts"],
        ends=lambda wc:    _DR_GENE_PARAMS[wc.drug]["ends"],
        strands=lambda wc: _DR_GENE_PARAMS[wc.drug]["strands"],
        results_dir=DR_DIR,
    shell:
        r"""
        set -euo pipefail
        Rscript scripts/dr_mutation_selection/02-make_initial_list.R \
            {params.drug} {params.genes} {params.starts} {params.ends} \
            {params.strands} {params.results_dir}
        """


# =============================================================================
# Step 11: DR – Apply Threshold × Promoter-Length Combination to Generate List1
# =============================================================================
# Applies convergence-event-number threshold and promoter-length filter to
# the annotated initial candidate list to produce list1.

rule dr_make_list1:
    input:
        initial=f"{DR_DIR}/{{drug}}/denovo_EvoResist_initial_list.txt",
    output:
        f"{DR_DIR}/{{drug}}/Threshold_{{threshold}}_Promoter_{{promoter}}_list1.tsv",
    params:
        drug="{drug}",
        genes=lambda wc:   _DR_GENE_PARAMS[wc.drug]["genes"],
        starts=lambda wc:  _DR_GENE_PARAMS[wc.drug]["starts"],
        ends=lambda wc:    _DR_GENE_PARAMS[wc.drug]["ends"],
        lof=lambda wc:     _DR_GENE_PARAMS[wc.drug]["lof"],
        strands=lambda wc: _DR_GENE_PARAMS[wc.drug]["strands"],
        results_dir=DR_DIR,
    shell:
        r"""
        set -euo pipefail
        Rscript scripts/dr_mutation_selection/04-make_thres_prom_combination.R \
            {params.drug} {params.genes} {params.starts} {params.ends} \
            {params.lof} {wildcards.threshold} {wildcards.promoter} \
            {params.strands} {params.results_dir}
        """


# =============================================================================
# Step 12: DR – Leave-One-Out Evaluation of a Variant List on Train or Test Split
# =============================================================================
# Computes per-isolate predictions, overall metrics (sensitivity, specificity,
# PPV, NPV with 95 % Wilson CIs), and per-variant leave-one-out deltas.
# Wildcard `listver` is `list1` or `list2`; `split` is `train` or `test`.

rule dr_loo_evaluate:
    input:
        variants=f"{DR_DIR}/{{drug}}/Threshold_{{threshold}}_Promoter_{{promoter}}_{{listver}}.tsv",
        id_list=lambda wc: (
            f"{DR_DIR}/{wc.drug}/id/train_70.txt"
            if wc.split == "train"
            else f"{DR_DIR}/{wc.drug}/id/test_30.txt"
        ),
    output:
        overall=f"{DR_DIR}/{{drug}}/Threshold_{{threshold}}_Promoter_{{promoter}}/{{split}}_{{listver}}/overall_metrics.tsv",
        per_var=f"{DR_DIR}/{{drug}}/Threshold_{{threshold}}_Promoter_{{promoter}}/{{split}}_{{listver}}/per_variant_analysis.tsv",
        preds=f"{DR_DIR}/{{drug}}/Threshold_{{threshold}}_Promoter_{{promoter}}/{{split}}_{{listver}}/isolate_predictions.tsv",
    params:
        snp_dir=_DR_SNP_ANNO_DIR,
        indel_dir=_DR_INDEL_ANNO_DIR,
        outdir=lambda wc: (
            f"{DR_DIR}/{wc.drug}/Threshold_{wc.threshold}_Promoter_{wc.promoter}"
            f"/{wc.split}_{wc.listver}"
        ),
    shell:
        r"""
        set -euo pipefail
        python3 scripts/dr_mutation_selection/03-leave_one_out.py \
            --variants_file {input.variants} \
            --id_list_file  {input.id_list} \
            --snp_dir       {params.snp_dir} \
            --indel_dir     {params.indel_dir} \
            --output_dir    {params.outdir}
        """


# =============================================================================
# Step 13: DR – Apply LOO Filtering Criteria to Generate Refined List2
# =============================================================================
# Reads per-variant leave-one-out results from the training-set evaluation of
# list1 and removes variants whose removal improves combined sensitivity +
# specificity, yielding a refined list2.

rule dr_make_list2:
    input:
        loo=f"{DR_DIR}/{{drug}}/Threshold_{{threshold}}_Promoter_{{promoter}}/train_list1/per_variant_analysis.tsv",
    output:
        list2=f"{DR_DIR}/{{drug}}/Threshold_{{threshold}}_Promoter_{{promoter}}_list2.tsv",
        removed=f"{DR_DIR}/{{drug}}/Threshold_{{threshold}}_Promoter_{{promoter}}_list1_removingrecords.tsv",
    params:
        drug="{drug}",
        genes=lambda wc:   _DR_GENE_PARAMS[wc.drug]["genes"],
        starts=lambda wc:  _DR_GENE_PARAMS[wc.drug]["starts"],
        ends=lambda wc:    _DR_GENE_PARAMS[wc.drug]["ends"],
        lof=lambda wc:     _DR_GENE_PARAMS[wc.drug]["lof"],
        strands=lambda wc: _DR_GENE_PARAMS[wc.drug]["strands"],
        results_dir=DR_DIR,
    shell:
        r"""
        set -euo pipefail
        Rscript scripts/dr_mutation_selection/05-loo_evaluate.R \
            {params.drug} {params.genes} {params.starts} {params.ends} \
            {params.lof} {wildcards.threshold} {wildcards.promoter} \
            {params.strands} {params.results_dir}
        """

