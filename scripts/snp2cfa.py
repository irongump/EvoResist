#!/usr/bin/env python3
"""
snp2cfa.py  --  Build a concatenated-FASTA alignment from per-sample files.

Usage:
    python snp2cfa.py <pos_file> <cfa_strain_file> [ref_fasta]

Arguments:
    pos_file         - Text file; one 1-based genomic position per line (sorted).
    cfa_strain_file  - Text file; one file path per line.
                       Each file can be:
                         .cfa / .fa / .fasta : full-genome FASTA consensus
                         .snp               : 3-column (pos\\tref\\talt)
                         other              : treated as SNP format
    ref_fasta        - Optional.  Path to the reference/ancestor FASTA.
                       If omitted the script looks for the last entry in
                       cfa_strain_file whose name starts with 'tb.ancestor'.

Outputs (derived from pos_file path  <dir>/<lineage>.pos):
    <dir>/<lineage>.fa                        - multi-FASTA alignment
    <dir>/../lineage_cfa/<lineage>_delete.pos - ordered position list used
                                                (for nodes_base_locus_iqtree.pl)
"""

import os
import sys


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def load_reference(fasta_path):
    """Return genome sequence (str, 0-indexed, uppercase) from a FASTA file."""
    seq_parts = []
    with open(fasta_path) as fh:
        for line in fh:
            if line.startswith(">"):
                continue
            seq_parts.append(line.strip().upper())
    return "".join(seq_parts)


def load_cfa_fasta(fasta_path):
    """Return {seq_name: sequence_str} from a (potentially multi-line) FASTA."""
    seqs = {}
    name = None
    buf = []
    with open(fasta_path) as fh:
        for line in fh:
            line = line.rstrip()
            if line.startswith(">"):
                if name is not None:
                    seqs[name] = "".join(buf)
                name = line[1:].strip().split()[0]
                buf = []
            else:
                buf.append(line)
    if name is not None:
        seqs[name] = "".join(buf)
    return seqs


def load_snp_file(snp_path):
    """Return {pos (int): alt (str)} mapping from a 3-column SNP file."""
    snps = {}
    with open(snp_path) as fh:
        for line in fh:
            line = line.strip()
            if not line:
                continue
            parts = line.split()
            if len(parts) >= 3:
                try:
                    pos = int(parts[0])
                    alt = parts[2].upper()
                    snps[pos] = alt
                except ValueError:
                    continue
    return snps


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main():
    if len(sys.argv) < 3:
        sys.exit(f"Usage: {sys.argv[0]} <pos_file> <cfa_strain_file> [ref_fasta]")

    pos_file = sys.argv[1]
    cfa_strain_file = sys.argv[2]
    provided_ref = sys.argv[3] if len(sys.argv) > 3 else None

    # -- Derive output paths from pos_file -----------------------------------
    pos_dir = os.path.dirname(os.path.abspath(pos_file))
    lineage = os.path.splitext(os.path.basename(pos_file))[0]  # strip .pos
    output_fa = os.path.join(pos_dir, f"{lineage}.fa")
    delete_pos_dir = os.path.normpath(os.path.join(pos_dir, "..", "lineage_cfa"))
    os.makedirs(delete_pos_dir, exist_ok=True)
    delete_pos_file = os.path.join(delete_pos_dir, f"{lineage}_delete.pos")

    # -- Read positions -------------------------------------------------------
    with open(pos_file) as fh:
        positions = [int(line.strip()) for line in fh if line.strip()]

    if not positions:
        sys.exit(f"Error: pos_file '{pos_file}' is empty.")

    # -- Read file list -------------------------------------------------------
    with open(cfa_strain_file) as fh:
        file_paths = [line.strip() for line in fh if line.strip()]

    if not file_paths:
        sys.exit(f"Error: cfa_strain_file '{cfa_strain_file}' is empty.")

    # -- Locate reference FASTA (for SNP-file fallback) ----------------------
    ref_genome = None
    if provided_ref:
        ref_genome = load_reference(provided_ref)
    else:
        # Last entry whose basename starts with 'tb.ancestor' is the ancestor
        for fp in reversed(file_paths):
            bn = os.path.basename(fp)
            if bn.startswith("tb.ancestor"):
                ref_genome = load_reference(fp)
                break
        if ref_genome is None:
            # Try default location
            default_ref = "data/tb.ancestor.fasta"
            if os.path.exists(default_ref):
                ref_genome = load_reference(default_ref)

    # -- Build alignment sequences -------------------------------------------
    records = []  # list of (name, sequence_str)

    for fp in file_paths:
        ext = os.path.splitext(fp)[1].lower()
        name = os.path.splitext(os.path.basename(fp))[0]

        if ext in (".cfa", ".fa", ".fasta"):
            # Full-genome FASTA: extract bases at each position (1-indexed)
            seqs = load_cfa_fasta(fp)
            # The FASTA may have one or more sequences; use the first (or match name)
            if not seqs:
                print(f"Warning: FASTA file '{fp}' is empty, skipping.", file=sys.stderr)
                continue
            seq = next(iter(seqs.values()))  # take first sequence
            seq_name = next(iter(seqs.keys()))
            bases = []
            for pos in positions:
                idx = pos - 1  # 0-indexed
                if idx < len(seq):
                    b = seq[idx].upper()
                    bases.append(b if b in "ACGT" else "N")
                else:
                    bases.append("N")
            records.append((seq_name, "".join(bases)))

        else:
            # SNP file (or unknown): need reference genome
            if ref_genome is None:
                print(
                    f"Warning: '{fp}' looks like a SNP file but no reference "
                    f"genome is available.  Skipping.",
                    file=sys.stderr,
                )
                continue
            if not os.path.exists(fp):
                print(f"Warning: '{fp}' not found, skipping.", file=sys.stderr)
                continue
            snps = load_snp_file(fp)
            bases = []
            for pos in positions:
                if pos in snps:
                    b = snps[pos].upper()
                    bases.append(b if b in "ACGT" else "N")
                else:
                    idx = pos - 1
                    if idx < len(ref_genome):
                        b = ref_genome[idx].upper()
                        bases.append(b if b in "ACGT" else "N")
                    else:
                        bases.append("N")
            records.append((name, "".join(bases)))

    if not records:
        sys.exit("Error: No valid sequences could be built for the alignment.")

    # -- Write alignment FASTA -----------------------------------------------
    with open(output_fa, "w") as out:
        for seq_name, seq in records:
            out.write(f">{seq_name}\n{seq}\n")

    print(
        f"snp2cfa: wrote {len(records)} sequences × {len(positions)} sites "
        f"→ {output_fa}",
        file=sys.stderr,
    )

    # -- Write delete.pos (ordered position list) ----------------------------
    with open(delete_pos_file, "w") as out:
        for pos in positions:
            out.write(f"{pos}\n")

    print(f"snp2cfa: wrote position list → {delete_pos_file}", file=sys.stderr)


if __name__ == "__main__":
    main()
