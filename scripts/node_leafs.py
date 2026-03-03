#!/usr/bin/env python3
"""
node_leafs.py  --  Annotate a branch-mutation file with per-node leaf counts.

Usage:
    python node_leafs.py <treefile> <db_mutation_file>

The db_mutation_file has tab-separated columns:
    node_id  col2  mutation  db  [...]

This script appends the leaf count (number of tip/leaf sequences that descend
from each node) as a new last column and prints the result to stdout.
"""

import sys
import re
from collections import defaultdict


# ---------------------------------------------------------------------------
# Newick parser (minimal)
# ---------------------------------------------------------------------------

def parse_newick(newick_str):
    """
    Parse a Newick string and return:
        parent_of  {node_name: parent_name or None}
        children   {node_name: [child_names]}
        is_leaf    {node_name: bool}
    """
    # Strip whitespace / trailing semicolon
    s = newick_str.strip().rstrip(";")

    parent_of = {}
    children = defaultdict(list)
    node_counter = [0]

    def _new_internal():
        node_counter[0] += 1
        return f"__internal_{node_counter[0]}"

    def _parse(pos, parent):
        """Recursive descent; returns the name of the parsed node."""
        if pos >= len(s):
            return None, pos

        if s[pos] == "(":
            # Internal node: parse children
            pos += 1  # skip '('
            node_name = None  # will be set after children
            child_names = []
            while pos < len(s) and s[pos] != ")":
                child_name, pos = _parse(pos, None)  # parent set below
                child_names.append(child_name)
                if pos < len(s) and s[pos] == ",":
                    pos += 1  # skip ','
            pos += 1  # skip ')'
            # Read optional node label
            label, pos = _read_label(pos)
            # Read optional branch length
            pos = _skip_branch_length(pos)
            node_name = label if label else _new_internal()
            for ch in child_names:
                children[node_name].append(ch)
                parent_of[ch] = node_name
            if node_name not in parent_of:
                parent_of[node_name] = parent
            return node_name, pos
        else:
            # Leaf node
            label, pos = _read_label(pos)
            pos = _skip_branch_length(pos)
            leaf_name = label if label else _new_internal()
            parent_of[leaf_name] = parent
            return leaf_name, pos

    def _read_label(pos):
        label = []
        while pos < len(s) and s[pos] not in (",", ")", ":", ";", "("):
            label.append(s[pos])
            pos += 1
        return "".join(label).strip(), pos

    def _skip_branch_length(pos):
        if pos < len(s) and s[pos] == ":":
            pos += 1
            while pos < len(s) and s[pos] not in (",", ")", ";", "("):
                pos += 1
        return pos

    root_name, _ = _parse(0, None)

    # Determine which nodes are leaves
    all_nodes = set(parent_of.keys()) | {root_name}
    is_leaf = {n: (n not in children or len(children[n]) == 0) for n in all_nodes}

    return parent_of, dict(children), is_leaf, root_name


# ---------------------------------------------------------------------------
# Count leaves under each node
# ---------------------------------------------------------------------------

def count_leaves(children, is_leaf, root):
    """Return {node_name: leaf_count} for all nodes."""
    leaf_count = {}

    def _count(node):
        if node in leaf_count:
            return leaf_count[node]
        if is_leaf.get(node, True):
            leaf_count[node] = 1
        else:
            total = 0
            for ch in children.get(node, []):
                total += _count(ch)
            leaf_count[node] = total
        return leaf_count[node]

    _count(root)
    # Ensure every leaf also has a count
    for node, flag in is_leaf.items():
        if flag and node not in leaf_count:
            leaf_count[node] = 1
    return leaf_count


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main():
    if len(sys.argv) < 3:
        sys.exit(f"Usage: {sys.argv[0]} <treefile> <db_mutation_file>")

    treefile = sys.argv[1]
    mutation_file = sys.argv[2]

    with open(treefile) as fh:
        newick = fh.read().strip()

    _, children, is_leaf, root = parse_newick(newick)
    leaf_count = count_leaves(children, is_leaf, root)

    with open(mutation_file) as fh:
        for line in fh:
            line = line.rstrip("\n")
            if not line:
                continue
            parts = line.split("\t")
            node = parts[0]
            count = leaf_count.get(node, 1)
            print(f"{line}\t{count}")


if __name__ == "__main__":
    main()
