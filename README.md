# EvoResist

Evolution-Guided Prioritization of Drug Resistance Mutations Enhances Molecular Prediction of Tuberculosis Drug Susceptibility.

---

## Pipeline Overview

The workflow consists of **19 steps** (Steps 1–7: convergent evolution analysis; Steps 8–19: DR mutation selection):

| Step | Rule name(s) | Description |
|------|-------------|-------------|
| 1 | `snp_calling` | SNP calling from FASTQ or SRA files (per sample) |
| 2 | `indel_calling` | INDEL calling from BAM files (per sample) |
| 3 | `build_tree` | Phylogenetic tree building with IQ-TREE (per lineage/sublineage) |
| 4 | `branch_mutations` | Branch mutations extraction (within lineage/sublineage)|
| 5 | `ancestor_mutations` | Ancestor mutations extraction (prior to lineage/sublineage diversification) |
| 6 | `merge_annotations` → `stat_convergent` → `filter_convergent` | Count convergent mutations by codon |
| 7 | `simulation` | GTR+Gamma simulation of mutations under a null distribution |
| 8 | `dr_train_test_split` | Stratify 70/30 train-test for sensitivity analysis to identify drug specific convergent threshold |
| 9 | `dr_filter_variants` | Include convergent SNPs or indels from drug resistance genes and promoter regions |
| 10 | `dr_initial_list` | Build annotated initial candidate variant list per drug |
| 11 | `dr_make_list1` | Apply threshold × promoter-length combination to generate list1 |
| 12 | `dr_loo_evaluate` | Leave-one-out evaluation of a variant list on train or test split |
| 13 | `dr_make_list2` | Apply LOO filtering criteria to generate refined list2 |
| 14 | `dr_select_best_threshold` | Select best convergence threshold per drug from threshold sweep (maximise MCC on test set) |
| 15 | `dr_select_best_promoter` | Select best promoter length per drug from promoter-length sweep (maximise MCC on test set) |
| 16 | `dr_final_evaluate` | Final evaluation on full cohort for EvoResist, WHO G1+G2, and WHO G1 catalogues |
| 17 | `dr_gain_evaluation` | Incremental gain of EvoResist variants relative to WHO G1+G2 baseline |
| 18 | `dr_lasso_analysis` | LASSO logistic regression for variant-phenotype association adjusted for lineage |
| 19 | `dr_compare_who_evoresist` | Compare EvoResist vs WHO G1 and G1+G2 performance (overall and stratified) |

---

## Installation

### 1. Clone the repository

```bash
git clone https://github.com/irongump/EvoResist.git
cd EvoResist
```

### 2. Create the conda environment

All required software is listed in `environment.yaml`.  
Install with [conda](https://docs.conda.io/) or [mamba](https://mamba.readthedocs.io/):

```bash
conda env create -f environment.yaml
conda activate evoresist
```

> **Note:** [mamba](https://mamba.readthedocs.io/) is strongly recommended as a faster drop-in replacement for conda:
> ```bash
> mamba env create -f environment.yaml
> ```

### 3. Install VarScan manually

VarScan is not available on conda.  Download the JAR from the [VarScan releases page](https://github.com/dkoboldt/varscan/releases) and place it at the path configured in `config/config.yaml` (default: `scripts/VarScan.v2.3.9.jar`):

```bash
wget -P scripts/ \
  https://github.com/dkoboldt/varscan/releases/download/2.3.9/VarScan.v2.3.9.jar
```

### 4. Install Picard manually

Picard is not available on conda.  Download the JAR from the [Picard releases page](https://github.com/broadinstitute/picard/releases) and place it at the path configured in `config/config.yaml` (default: `scripts/picard.jar`):

```bash
wget -P scripts/ \
  https://github.com/broadinstitute/picard/releases/latest/download/picard.jar
```

### 5. Install GATK3 manually

GATK3 (`GenomeAnalysisTK.jar`) requires a license and manual download from the [GATK archive](https://console.cloud.google.com/storage/browser/gatk-software/package-archive/gatk).  Place the JAR at the path configured in `config/config.yaml` (default: `scripts/GenomeAnalysisTK.jar`).  **Java 8 (JDK 1.8)** is required for GATK3 compatibility.

### Required software summary

| Tool | Version | Purpose |
|------|---------|---------|
| [Snakemake](https://snakemake.readthedocs.io/) | ≥ 7.0 | Workflow management |
| [BWA](http://bio-bwa.sourceforge.net/) | ≥ 0.7.17 | Short-read alignment |
| [SAMtools](http://www.htslib.org/) | ≥ 1.13 | SAM/BAM processing |
| [sickle](https://github.com/najoshi/sickle) | ≥ 1.33 | Quality trimming |
| [SRA Toolkit](https://github.com/ncbi/sra-tools) | ≥ 3.0 | SRA → FASTQ conversion |
| [IQ-TREE 2](http://www.iqtree.org/) | ≥ 2.2 | Phylogenetic tree inference |
| [GNU parallel](https://www.gnu.org/software/parallel/) | ≥ 20210722 | Parallel annotation |
| [bcftools](http://www.htslib.org/) | ≥ 1.13 | VCF normalization and filtering (INDEL calling) |
| [VarScan](https://github.com/dkoboldt/varscan) | 2.3.9 | Variant calling (manual install) |
| [Picard](https://broadinstitute.github.io/picard/) | ≥ 2.27 | Read group addition for INDEL calling (manual install) |
| [GATK3](https://gatk.broadinstitute.org/) | 3.x | Local realignment and INDEL calling (manual install) |
| [Java](https://openjdk.org/) | 8 (GATK3) / ≥ 11 (others) | Required by VarScan and GATK3 |
| [Perl](https://www.perl.org/) | ≥ 5.26 | Annotation scripts |
| [R](https://www.r-project.org/) | ≥ 4.0 | Statistical analysis |
| [Python](https://www.python.org/) | ≥ 3.8 | Simulation and filtering |
| numpy, pandas, biopython, tqdm | latest | Python dependencies |
| r-dplyr, r-ggplot2, r-ggpubr, r-data.table, … | latest | R dependencies |

---

## Configuration

Edit `config/config.yaml` before running.  Key options:

### Steps 1–7 (convergent evolution analysis)

| Option | Default | Description |
|--------|---------|-------------|
| `outdir` | `"output"` | Base directory for all pipeline outputs |
| `stop_at` | `"all"` | Stop the pipeline at a specific step (see table above for valid step names) |
| `input_type` | `"auto"` | Input data format: `auto`, `fastq`, or `sra` |
| `fastq_dir` | *(none)* | **Required when `input_type: fastq`** — path to the directory containing pre-existing per-sample FASTQ files. Ignored for `sra` and `auto` modes. |
| `lineage` | *(none — all)* | Process only the named lineage or list of lineages (e.g. `"Lineage1.1.A"` or `["Lineage1.1.A","Lineage2.3.4"]`). Each must correspond to a `<strain_ids_dir>/<lineage>_strain.txt` file. When omitted, all lineages are processed. |

### Steps 8–9 (DR mutation selection)

| Option | Default | Description |
|--------|---------|-------------|
| `dr_results_dir` | `"dr_results"` | Base directory for DR mutation selection outputs |
| `snp_anno_dir` | *(required)* | Directory containing per-sample SNP annotation files named `{sample}.ano` (tab-separated, columns: position, ref, alt). **Must be set** when running steps 8–19. |
| `indel_anno_dir` | *(required)* | Directory containing per-sample indel annotation files named `{sample}.indel.ano` (same format). **Must be set** when running steps 8–19. |
| `all_indel_file` | `"data/all_indel_100k.txt.gz"` | Cohort-wide indel file used to build per-drug indel candidate lists. Accepts `.gz` compressed or plain text. |
| `dr_convergent_snp_file` | *(step 6 output)* | Convergent SNP file for DR analysis. Defaults to `<outdir>/lineage_ann/all_ann_convergent_flt.txt` (step 6 output). Override to use a pre-existing file (e.g. `data/all_ann_convergent_flt.txt.gz`). |
| `dr_best_promoters` | *(500 for all drugs)* | Per-drug best promoter length (bp) for final evaluation and LASSO steps. Determined after running step 11 (`dr_best_promoter`). Override per drug as a YAML dict (e.g. `{RIF: 400, INH: 300}`). |

`stop_at` accepts either the step number (`step1`–`step14`/`dr_prep`/`dr_selection`/…) or the rule name.  
Valid `stop_at` values:

| Value | Alias | What is produced |
|-------|-------|-----------------|
| `step1` | `snp_calling` | Per-sample SNP files |
| `step2` | `indel_calling` | Per-sample INDEL VCF files |
| `step3` | `build_tree` | Per-lineage phylogenetic trees |
| `step4` | `branch_mutations` | Per-lineage branch mutation annotations |
| `step5` | `ancestor_mutations` | Ancestor mutation annotation |
| `step6` | `merge_annotations` | Merged and filtered convergent mutation list |
| `step7` | `simulation` / `all` | GTR simulation null distribution |
| `step8` | `dr_prep` | DR data preparation (train/test split, filtered candidates, initial lists) |
| `step9` | `dr_selection` | Full DR sweep: list1→LOO→list2→train/test evaluation for all threshold × promoter combinations |
| `step10` | `dr_best_threshold` | Best-threshold selection files (`best_thresholds.tsv`, `best_thresholds.sh`) |
| `step11` | `dr_best_promoter` | Best-promoter selection files (`best_promoters.tsv`, `best_promoters.sh`) |
| `step12` | `dr_final_evaluate` | Final full-cohort evaluation for EvoResist, WHO G1+G2, and WHO G1 |
| `step13` | `dr_lasso` | Incremental-gain ranking + LASSO logistic regression outputs |
| `step14` | `dr_compare` | WHO vs EvoResist comparison tables (overall + stratified) |

`input_type` values:
- `auto` – look for FASTQ files in `<outdir>/fastq/`; if absent, convert from SRA into the same directory
- `fastq` – FASTQ files already exist at the path given by `fastq_dir` (required); only samples whose FASTQ files are present in `fastq_dir` are processed (others are skipped with a warning)
- `sra` – input is SRA format; convert to FASTQ into `<outdir>/fastq/` before processing

Options can also be overridden on the command line with `--config`:

```bash
# Run only SNP calling (step 1):
snakemake --cores 8 --configfile config/config.yaml --config stop_at=step1

# Use a custom output directory:
snakemake --cores 8 --configfile config/config.yaml --config outdir=results

# Start from SRA files:
snakemake --cores 8 --configfile config/config.yaml --config input_type=sra

# Start from existing FASTQ files in a custom directory:
snakemake --cores 8 --configfile config/config.yaml \
    --config input_type=fastq fastq_dir=/path/to/my/fastqs

# Process only Lineage1.1.A from existing FASTQ files:
snakemake --cores 8 --configfile config/config.yaml \
    --config input_type=fastq fastq_dir=/path/to/my/fastqs lineage=Lineage1.1.A

# Process Lineage1.1.A and Lineage2.3.4 (edit config.yaml for multi-value):
#   lineage: ["Lineage1.1.A", "Lineage2.3.4"]
```

---

## Usage

```bash
# Run the full convergent evolution analysis (steps 1–7):
snakemake --cores <N> --configfile config/config.yaml

# Dry-run to check which jobs will be executed:
snakemake --cores <N> --configfile config/config.yaml -n

# Run the DR mutation selection analysis (steps 8–9).
# Requires snp_anno_dir and indel_anno_dir to be configured.
# If the step 6 output already exists, it is used directly without re-running
# steps 1–7. To use a pre-existing convergent SNP file instead, set
# dr_convergent_snp_file in config.yaml.
snakemake --cores <N> --configfile config/config.yaml --config stop_at=dr_selection

# Submit jobs to an HPC cluster (SLURM example):
snakemake --cores <N> --configfile config/config.yaml \
    --cluster "sbatch --mem={resources.mem_mb}M --cpus-per-task={threads}"
```

### DR analysis – required data files

Before running steps 8–19, ensure the following files are present:

| File | Required from step | Description |
|------|-------------------|-------------|
| `data/{drug}_sample_list.txt` | 8 | Per-drug sample list with columns `Run` (sample ID) and `pheno` (R/S phenotype). One file per drug (RIF, INH, EMB, PZA, LFX, MFX, BDQ, AMK, STM, ETO, KAN, CAP, LZD). |
| `{dr_results_dir}/{drug}/WHO_list/WHO_list_allGroup.txt` | 10 | WHO variant catalogue for each drug (headerless TSV). Must be placed in the correct subdirectory before running `dr_initial_list`. |
| `{dr_results_dir}/{drug}/WHO_list/WHO_G1_2_withcolnames.txt` | 16 | WHO Grade 1+2 variant list with column headers. Required for `dr_final_evaluate`, `dr_gain_evaluation`, and `dr_compare_who_evoresist`. |
| `{dr_results_dir}/{drug}/WHO_list/WHO_G1_withcolnames.txt` | 16 | WHO Grade 1 variant list with column headers. Required for `dr_final_evaluate` and `dr_compare_who_evoresist`. |
| `data/all_indel_100k.txt.gz` | 9 | Cohort-wide indel file (configurable via `all_indel_file`). |
| Per-sample SNP files | 8 | Files named `{sample}.ano` in the directory set by `snp_anno_dir`. |
| Per-sample indel files | 8 | Files named `{sample}.indel.ano` in the directory set by `indel_anno_dir`. |

### Steps 14–19 – automated final analysis

After the full DR sweep (step 9 / `dr_selection`), the remaining steps run automatically via Snakemake.
Note: the pipeline overview table above numbers all Snakemake rules as steps 1–19; the `stop_at`
values below (`step10`–`step14`) group several of those rules into logical checkpoints.

**Select best threshold** (`stop_at=dr_best_threshold` / `step10`):

```bash
snakemake --cores <N> --configfile config/config.yaml --config stop_at=dr_best_threshold
# Produces: dr_results/best_thresholds.tsv  and  dr_results/best_thresholds.sh
```

The best thresholds are already hardcoded in `Snakefile` (`_DR_BEST_THRESHOLDS`). The rule writes the selection report but the hardcoded values drive the promoter sweep in step 9.

**Select best promoter** (`stop_at=dr_best_promoter` / `step11`):

```bash
snakemake --cores <N> --configfile config/config.yaml --config stop_at=dr_best_promoter
# Produces: dr_results/best_promoters.tsv  and  dr_results/best_promoters.sh
```

After this step, update `dr_best_promoters` in `config/config.yaml` with the values from `best_promoters.tsv` before running the final analysis steps.

**Final evaluation, gain, LASSO, comparison** (`stop_at=dr_compare` / `step14`):

```bash
snakemake --cores <N> --configfile config/config.yaml --config stop_at=dr_compare
```

This requires `dr_best_promoters` to be configured and the WHO catalogue files (`WHO_G1_2_withcolnames.txt`, `WHO_G1_withcolnames.txt`) to be present in each `{dr_results_dir}/{drug}/WHO_list/` directory.

---

## Output structure

```
<outdir>/
├── fastq/          # FASTQ files (after SRA conversion, if applicable)
├── bam/            # Aligned and sorted BAM files
├── indel/          # Per-sample INDEL VCF files (raw, normalized, and pass-filtered)
├── snv/            # Per-sample SNP / pileup files
├── cfa/            # CFA consensus files
├── forup/          # forup files
├── lineage_tree/   # IQ-TREE output (tree, state files)
├── lineage_ann/    # Per-lineage and merged annotation files
└── simulation/     # GTR simulation outputs (CSV)

<dr_results_dir>/           # DR mutation selection outputs (default: dr_results/)
└── {drug}/                 # One directory per drug
    ├── id/
    │   ├── train_70.txt    # Training-set sample IDs (70 %)
    │   └── test_30.txt     # Test-set sample IDs (30 %)
    ├── WHO_list/
    │   └── WHO_list_allGroup.txt   # User-provided WHO catalogue (input)
    ├── denovo_snp_2.txt            # Drug-region-filtered convergent SNPs
    ├── denovo_indel_2.txt          # Drug-region-filtered indels
    ├── denovo_EvoResist_initial_list.txt   # Annotated initial candidate list
    ├── Threshold_{T}_Promoter_{P}_list1.tsv        # list1 (threshold × promoter)
    ├── Threshold_{T}_Promoter_{P}_list2.tsv        # list2 (after LOO filtering)
    ├── Threshold_{T}_Promoter_{P}_list1_removingrecords.tsv
    └── Threshold_{T}_Promoter_{P}/
        ├── train_list1/    # LOO evaluation of list1 on training set
        │   ├── overall_metrics.tsv
        │   ├── per_variant_analysis.tsv
        │   └── isolate_predictions.tsv
        ├── train_list2/    # LOO evaluation of list2 on training set
        └── test_list2/     # Evaluation of list2 on test set
├── best_thresholds.tsv         # Best convergence threshold per drug (step 14)
├── best_thresholds.sh          # Bash-sourceable threshold variables
├── best_promoters.tsv          # Best promoter length per drug (step 15)
├── best_promoters.sh           # Bash-sourceable promoter variables
└── {drug}/
    ├── Final_Evaluate/
    │   ├── EvoResist/          # Full-cohort evaluation of EvoResist final list (step 16)
    │   │   ├── overall_metrics.tsv
    │   │   ├── per_variant_analysis.tsv
    │   │   └── isolate_predictions.tsv
    │   ├── G1_2/               # Full-cohort evaluation of WHO G1+G2 (step 16)
    │   └── G1/                 # Full-cohort evaluation of WHO G1 (step 16)
    ├── gain_EvoResist.txt      # Incremental gain ranking vs WHO G1+G2 (step 17)
    └── lasso/                  # LASSO logistic regression outputs (step 18)
        ├── data_summary.tsv
        ├── cv_performance.tsv
        ├── model_coefficients.tsv
        └── refit_coefficients.tsv
comparison/                     # WHO vs EvoResist comparison (step 19)
├── overall/
│   └── comparison_formatted.tsv
├── by_lineage/
│   ├── comparison_formatted.tsv
│   └── comparison_stratified_formatted.tsv
└── by_pdst_method/
    ├── comparison_formatted.tsv
    └── comparison_stratified_formatted.tsv
```
