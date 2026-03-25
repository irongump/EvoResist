# EvoResist

Leveraging Convergent Evolution to Prioritize Antibiotic Resistance Mutations in *Mycobacterium tuberculosis*.

---

## Pipeline Overview

The workflow consists of **6 steps**:

| Step | Rule name | Description |
|------|-----------|-------------|
| 1 | `snp_calling` | SNP calling from FASTQ or SRA files (per sample) |
| 2 | `build_tree` | Phylogenetic tree building with IQ-TREE (per lineage) |
| 3 | `branch_mutations` | Branch mutation extraction from tree nodes (per lineage) |
| 4 | `ancestor_mutations` | Ancestor mutation extraction (global) |
| 5 | `merge_annotations` → `stat_convergent` → `filter_convergent` | Merge annotations and count convergent mutations (global) |
| 6 | `simulation` | GTR+Gamma simulation for null distribution (global) |

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
| [VarScan](https://github.com/dkoboldt/varscan) | 2.3.9 | Variant calling (manual install) |
| [Java](https://openjdk.org/) | ≥ 11 | Required by VarScan |
| [Perl](https://www.perl.org/) | ≥ 5.26 | Annotation scripts |
| [R](https://www.r-project.org/) | ≥ 4.0 | Statistical analysis |
| [Python](https://www.python.org/) | ≥ 3.8 | Simulation and filtering |
| numpy, pandas, biopython, tqdm | latest | Python dependencies |
| r-dplyr, r-ggplot2, r-ggpubr, r-data.table, … | latest | R dependencies |

---

## Configuration

Edit `config/config.yaml` before running.  Key options:

| Option | Default | Description |
|--------|---------|-------------|
| `outdir` | `"output"` | Base directory for all pipeline outputs |
| `stop_at` | `"all"` | Stop the pipeline at a specific step (see table above for valid step names) |
| `input_type` | `"auto"` | Input data format: `auto`, `fastq`, or `sra` |
| `fastq_dir` | *(none)* | **Required when `input_type: fastq`** — path to the directory containing pre-existing per-sample FASTQ files. Ignored for `sra` and `auto` modes. |

`stop_at` accepts either the step number (`step1`–`step6`) or the rule name.  
`input_type` values:
- `auto` – look for FASTQ files in `<outdir>/fastq/`; if absent, convert from SRA into the same directory
- `fastq` – FASTQ files already exist at the path given by `fastq_dir` (required); SRA download is skipped
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
```

---

## Usage

```bash
# Run the full pipeline (6 steps):
snakemake --cores <N> --configfile config/config.yaml

# Dry-run to check which jobs will be executed:
snakemake --cores <N> --configfile config/config.yaml -n

# Submit jobs to an HPC cluster (SLURM example):
snakemake --cores <N> --configfile config/config.yaml \
    --cluster "sbatch --mem={resources.mem_mb}M --cpus-per-task={threads}"
```

---

## Output structure

```
<outdir>/
├── fastq/          # FASTQ files (after SRA conversion, if applicable)
├── bam/            # Aligned and sorted BAM files
├── snv/            # Per-sample SNP / pileup files
├── cfa/            # CFA consensus files
├── forup/          # forup files
├── lineage_tree/   # IQ-TREE output (tree, state files)
├── lineage_ann/    # Per-lineage and merged annotation files
└── simulation/     # GTR simulation outputs (CSV)
```
