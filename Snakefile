"""
EvoResist Snakemake Pipeline (v1.0)
====================================
Leveraging Convergent Evolution to Prioritize Antibiotic Resistance
Mutations in Mycobacterium tuberculosis.

Pipeline Steps:
  1. SNP calling from FASTQ files (per sample)
  2. Phylogenetic tree building with IQ-TREE (per lineage)
  3. Branch mutation extraction from tree nodes (per lineage)
  4. Ancestor mutation extraction (global)
  5. Merge annotations and count convergent mutations (global)
  6. GTR simulation for null distribution (global)
  7. FDR analysis and visualization (global)

Usage:
  snakemake --cores <N> --configfile config/config.yaml
  snakemake --cores <N> --configfile config/config.yaml --cluster "sbatch ..."
"""

import os

configfile: "config/config.yaml"

OUTPUT_DIR = config.get("output_dir", "output")
if os.sep in OUTPUT_DIR or "/" in OUTPUT_DIR:
    raise ValueError(
        f"output_dir must be a single-level directory name (e.g. 'output'), "
        f"got: '{OUTPUT_DIR}'"
    )

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
# Final target
# =============================================================================
rule all:
    input:
        f"{OUTPUT_DIR}/simulation/final_convergent_mutations_statistics.csv",
        f"{OUTPUT_DIR}/simulation/Figure_Null_Distribution_Barplot.pdf",

# =============================================================================
# Step 1: Build Phylogenetic Tree (per lineage)
# =============================================================================
# using different depth parameter to filter samples (e.g. 30x, 50x, 80x)
# Collects SNP positions across all samples in a lineage, generates a
# concatenated alignment, and builds a maximum-likelihood tree with IQ-TREE.

rule build_tree:
    input:
        snps=lambda wc: expand(OUTPUT_DIR + "/snv/{sample}.snp",
                               sample=get_samples(wc.lineage)),
        cfas=lambda wc: expand(OUTPUT_DIR + "/cfa/{sample}.cfa",
                               sample=get_samples(wc.lineage)),
        strain_list=os.path.join(config["strain_ids_dir"], "{lineage}.txt"),
        anc_cfa=config["ancestor_concat_fasta"],
        ref=config["reference"],
    output:
        tree=f"{OUTPUT_DIR}/lineage_tree/{{lineage}}_btp.treefile",
        state=f"{OUTPUT_DIR}/lineage_tree/{{lineage}}_btp.state",
    params:
        lineage="{lineage}",
        tree_dir=f"{OUTPUT_DIR}/lineage_tree",
        cfa_dir=f"{OUTPUT_DIR}/cfa",
        snv_dir=f"{OUTPUT_DIR}/snv",
        depth_thres=config["depth_thres"],
        depth_file=config.get("depth_file", ""),
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

        flt_strain_file="${{tree_dir}}/${{lineage}}_depth_{params.depth_thres}.txt"
        Rscript scripts/sample_flt.R {input.strain_list} {params.depth_file} {params.depth_thres}
        mv "$(basename {input.strain_list} .txt)_depth_{params.depth_thres}.txt" "$flt_strain_file"

        while IFS= read -r strain; do
            snp_file="${{snv_dir}}/${{strain}}.snp"
            if [ -f "$snp_file" ]; then
                awk '{{print $1}}' "$snp_file" >> "$pos_file"
                echo "$strain" >> "$snp_strain_file"
            else
                echo "Warning: ${{strain}} snp file not found, skipping." >&2
            fi
        done < "$flt_strain_file"

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
        tree=f"{OUTPUT_DIR}/lineage_tree/{{lineage}}_btp.treefile",
        state=f"{OUTPUT_DIR}/lineage_tree/{{lineage}}_btp.state",
        low_ebr=config["low_ebr_file"],
        ref=config["reference"],
    output:
        ann=f"{OUTPUT_DIR}/lineage_ann/{{lineage}}.ann",
    params:
        lineage="{lineage}",
        tree_dir=f"{OUTPUT_DIR}/lineage_tree",
        ann_dir=f"{OUTPUT_DIR}/lineage_ann",
        cfa_dir=f"{OUTPUT_DIR}/lineage_cfa",
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
        ann=f"{OUTPUT_DIR}/lineage_ann/ancestor.ann",
    params:
        ann_dir=f"{OUTPUT_DIR}/lineage_ann",
        anc_tmp=f"{OUTPUT_DIR}/ancestor_tmp",
    threads: config["threads_annotation"]
    resources:
        mem_mb=20000,
    shell:
        r"""
        set -euo pipefail
        ann_dir={params.ann_dir}
        low_ebr={input.low_ebr}
        mkdir -p "$ann_dir"

        # Generate per-node SNP files for ancestor nodes
        python scripts/getrefbase_per_node.py \
            {input.mutation_file} {params.anc_tmp} {input.ref}

        # Annotate each node's mutations in parallel
        find "{params.anc_tmp}" -maxdepth 1 -name "*snp" -print0 | \
            parallel -0 -j {threads} '
                python scripts/remove_low_ebr.py "'"$low_ebr"'" "{{}}" \
                    > "{{= s/\.snp$// =}}_rle.snp" &&
                perl scripts/mtbc_translate/0_MTBC_Annotation_mtbc_4411532_corrected.pl \
                    "{{= s/\.snp$// =}}_rle.snp" \
                    > "{{= s/\.snp$// =}}.ann" &&
                sed -i "/^$/d" "{{= s/\.snp$// =}}.ann"
            '

        cat {params.anc_tmp}/*.ann > "$ann_dir/ancestor.ann"

        # Clean up
        rm -rf {params.anc_tmp}
        """


# =============================================================================
# Step 5: Merge Annotations and Count Convergent Mutations
# =============================================================================
# Merges all lineage and ancestor annotations, then counts convergent events.

rule merge_annotations:
    input:
        lineage_anns=expand(f"{OUTPUT_DIR}/lineage_ann/{{lineage}}.ann",
                            lineage=LINEAGES),
        ancestor_ann=f"{OUTPUT_DIR}/lineage_ann/ancestor.ann",
    output:
        all_ann=f"{OUTPUT_DIR}/lineage_ann/all_ann.txt",
    shell:
        r"""
        cat {input.lineage_anns} > {output.all_ann}
        cat {input.ancestor_ann} >> {output.all_ann}
        """


rule stat_convergent:
    input:
        all_ann=f"{OUTPUT_DIR}/lineage_ann/all_ann.txt",
    output:
        convergent=f"{OUTPUT_DIR}/lineage_ann/all_ann_convergent.txt",
    params:
        ann_dir=f"{OUTPUT_DIR}/lineage_ann",
    shell:
        r"""
        cd {params.ann_dir}
        Rscript ../../scripts/stat_convergent.R
        """


rule filter_convergent:
    input:
        convergent=f"{OUTPUT_DIR}/lineage_ann/all_ann_convergent.txt",
        snp_freq=config["snp_freq_file"],
        repeat_region=config["repeat_region_file"],
        mobile_element=config["mobile_element_file"],
    output:
        filtered=f"{OUTPUT_DIR}/lineage_ann/all_ann_convergent_flt.txt",
    params:
        ann_dir=f"{OUTPUT_DIR}/lineage_ann",
    shell:
        r"""
        cd {params.ann_dir}
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
        filtered=f"{OUTPUT_DIR}/lineage_ann/all_ann_convergent_flt.txt",
        ref=config["reference"],
    output:
        raw_sim=f"{OUTPUT_DIR}/simulation/simulated_mutations_raw_GTR_Gamma.csv",
        null_dist=f"{OUTPUT_DIR}/simulation/null_mutation_df_GTR.csv",
        expected=f"{OUTPUT_DIR}/simulation/expected_null_distribution_GTR_Gamma.csv",
    params:
        sim_dir=f"{OUTPUT_DIR}/simulation",
    shell:
        r"""
        mkdir -p {params.sim_dir}
        cd {params.sim_dir}
        python ../../scripts/simulation_GTR_gamma.py
        """


# =============================================================================
# Step 7: FDR Analysis and Visualization
# =============================================================================
# Calculates empirical p-values and FDR for convergent mutations, and
# generates publication-ready plots.

rule fdr_analysis:
    input:
        filtered=f"{OUTPUT_DIR}/lineage_ann/all_ann_convergent_flt.txt",
        raw_sim=f"{OUTPUT_DIR}/simulation/simulated_mutations_raw_GTR_Gamma.csv",
    output:
        stats=f"{OUTPUT_DIR}/simulation/final_convergent_mutations_statistics.csv",
    params:
        sim_dir=f"{OUTPUT_DIR}/simulation",
    shell:
        r"""
        cd {params.sim_dir}
        Rscript ../../scripts/fdr4simulation.R
        """


rule plot_results:
    input:
        expected=f"{OUTPUT_DIR}/simulation/expected_null_distribution_GTR_Gamma.csv",
    output:
        pdf=f"{OUTPUT_DIR}/simulation/Figure_Null_Distribution_Barplot.pdf",
        png=f"{OUTPUT_DIR}/simulation/Figure_Null_Distribution_Barplot.png",
    params:
        sim_dir=f"{OUTPUT_DIR}/simulation",
    shell:
        r"""
        cd {params.sim_dir}
        Rscript ../../scripts/plot_GTR_simulation.R
        """
