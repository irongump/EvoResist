"""
EvoResist Snakemake Pipeline (v2.0)
====================================
Leveraging Convergent Evolution to Prioritize Antibiotic Resistance
Mutations in Mycobacterium tuberculosis.

Pipeline Steps:
  1. SNP calling from FASTQ or SRA files (per sample)
  2. Phylogenetic tree building with IQ-TREE (per lineage)
  3. Branch mutation extraction from tree nodes (per lineage)
  4. Ancestor mutation extraction (global)
  5. Merge annotations and count convergent mutations (global)
  6. GTR simulation for null distribution (global)

Key configuration options (config/config.yaml):
  outdir      - Base output directory (default: "output").
                Recommended: single-level relative path so that shell scripts
                using ../../ back-references resolve correctly.
  stop_at     - Stop after a specific step; valid values:
                  step1/snp_calling, step2/build_tree,
                  step3/branch_mutations, step4/ancestor_mutations,
                  step5/merge_annotations, step6/simulation/all
  input_type  - Input data format: auto (default), fastq, or sra.
                When set to "fastq", fastq_dir must also be set.
  fastq_dir   - Directory with pre-existing per-sample FASTQ files.
                Required when input_type is "fastq"; ignored otherwise.

Usage:
  snakemake --cores <N> --configfile config/config.yaml
  snakemake --cores <N> --configfile config/config.yaml --cluster "sbatch ..."
  # Run only up to SNP calling:
  snakemake --cores <N> --configfile config/config.yaml --config stop_at=step1
  # Use a custom output directory:
  snakemake --cores <N> --configfile config/config.yaml --config outdir=results
"""

import os

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
LINEAGES, = glob_wildcards(os.path.join(config["strain_ids_dir"], "{lineage}.txt"))
LINEAGES = sorted(LINEAGES)

def get_samples(lineage):
    """Read sample IDs from a lineage strain list file."""
    filepath = os.path.join(config["strain_ids_dir"], f"{lineage}.txt")
    with open(filepath) as f:
        return [line.strip() for line in f if line.strip()]

ALL_SAMPLES = []
for _lin in LINEAGES:
    ALL_SAMPLES.extend(get_samples(_lin))
ALL_SAMPLES = sorted(set(ALL_SAMPLES))

# =============================================================================
# Map each stop_at value to its corresponding output files
# =============================================================================
_SIMULATION_OUTPUTS = [
    f"{OUTDIR}/simulation/simulated_mutations_raw_GTR_Gamma.csv",
    f"{OUTDIR}/simulation/null_mutation_df_GTR.csv",
    f"{OUTDIR}/simulation/expected_null_distribution_GTR_Gamma.csv",
]

_STOP_STEP_MAP = {
    "step1":              expand(f"{OUTDIR}/snv/{{sample}}.snp", sample=ALL_SAMPLES),
    "snp_calling":        expand(f"{OUTDIR}/snv/{{sample}}.snp", sample=ALL_SAMPLES),
    "step2":              expand(f"{OUTDIR}/lineage_tree/{{lineage}}_btp.treefile", lineage=LINEAGES),
    "build_tree":         expand(f"{OUTDIR}/lineage_tree/{{lineage}}_btp.treefile", lineage=LINEAGES),
    "step3":              expand(f"{OUTDIR}/lineage_ann/{{lineage}}.ann", lineage=LINEAGES),
    "branch_mutations":   expand(f"{OUTDIR}/lineage_ann/{{lineage}}.ann", lineage=LINEAGES),
    "step4":              [f"{OUTDIR}/lineage_ann/ancestor.ann"],
    "ancestor_mutations": [f"{OUTDIR}/lineage_ann/ancestor.ann"],
    "step5":              [f"{OUTDIR}/lineage_ann/all_ann_convergent_flt.txt"],
    "merge_annotations":  [f"{OUTDIR}/lineage_ann/all_ann_convergent_flt.txt"],
    "step6":              _SIMULATION_OUTPUTS,
    "simulation":         _SIMULATION_OUTPUTS,
    "all":                _SIMULATION_OUTPUTS,
}

if STOP_AT not in _STOP_STEP_MAP:
    raise ValueError(
        f"Invalid stop_at value '{STOP_AT}'. "
        f"Valid options: {sorted(set(_STOP_STEP_MAP.keys()))}"
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
# Step 2: Build Phylogenetic Tree (per lineage)
# =============================================================================
# Collects SNP positions across all samples in a lineage, generates a
# concatenated alignment, and builds a maximum-likelihood tree with IQ-TREE.

rule build_tree:
    input:
        snps=lambda wc: expand(f"{OUTDIR}/snv/{{sample}}.snp",
                               sample=get_samples(wc.lineage)),
        cfas=lambda wc: expand(f"{OUTDIR}/cfa/{{sample}}.cfa",
                               sample=get_samples(wc.lineage)),
        strain_list=os.path.join(config["strain_ids_dir"], "{lineage}.txt"),
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
# Step 3: Branch Mutation Extraction (per lineage)
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
# Step 4: Ancestor Mutation Extraction (global)
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
# Step 5: Merge Annotations and Count Convergent Mutations
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
# Step 6: GTR Simulation for Null Distribution
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
