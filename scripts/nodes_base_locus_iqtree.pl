#!/usr/bin/perl
# nodes_base_locus_iqtree.pl
#
# Extract per-branch SNP mutations from IQ-TREE's marginal ancestral state
# reconstruction and write them in the standard db format.
#
# Usage (run from the lineage_tree directory):
#   perl nodes_base_locus_iqtree.pl \
#       <treefile>       \
#       <delete_pos>     \   # ordered list of genomic positions (1-based)
#       <state_file>     \   # IQ-TREE .state file
#       <alignment_fa>   \   # alignment FASTA used for tree building
#       <output_db>      \   # output: branch mutations in db format
#       <output_homo>        # output: homoplasy mutations (repeated across branches)
#
# Output (.db) format per line:
#   node_id TAB depth TAB POS_ALT TAB db
#
# where POS is the 1-based genomic position (from delete_pos) and ALT is the
# derived base at that node.

use strict;
use warnings;

# --------------------------------------------------------------------------
# Arguments
# --------------------------------------------------------------------------
my ($treefile, $delete_pos_file, $state_file, $aln_fa, $out_db, $out_homo) = @ARGV;

die "Usage: $0 treefile delete_pos state_file aln_fa out_db out_homo\n"
    unless defined $out_homo;

# --------------------------------------------------------------------------
# Read genomic positions (delete_pos = ordered list of positions used)
# --------------------------------------------------------------------------
my @positions;
if (-e $delete_pos_file) {
    open my $fh, '<', $delete_pos_file or die "Cannot open $delete_pos_file: $!";
    while (<$fh>) {
        chomp;
        push @positions, $_ if /^\d+$/;
    }
    close $fh;
}

# --------------------------------------------------------------------------
# Read alignment FASTA → leaf states
# --------------------------------------------------------------------------
my %leaf_seq;   # {sample_name => sequence_str}
{
    open my $fh, '<', $aln_fa or die "Cannot open $aln_fa: $!";
    my ($name, @parts);
    while (<$fh>) {
        chomp;
        if (/^>(\S+)/) {
            if (defined $name) { $leaf_seq{$name} = join('', @parts); }
            $name = $1;  @parts = ();
        } else {
            push @parts, $_;
        }
    }
    if (defined $name) { $leaf_seq{$name} = join('', @parts); }
    close $fh;
}

my $n_sites = @positions ? scalar(@positions) : do {
    my ($first_seq) = values %leaf_seq;
    defined $first_seq ? length($first_seq) : 0;
};

# If no positions file, create fake 1..n_sites positions
unless (@positions) {
    @positions = (1..$n_sites);
}

# --------------------------------------------------------------------------
# Parse Newick tree → parent_of, children, depth_from_root
# --------------------------------------------------------------------------
my %parent_of;
my %children;

{
    open my $fh, '<', $treefile or die "Cannot open $treefile: $!";
    my $nwk = do { local $/; <$fh> };
    close $fh;
    $nwk =~ s/\s+//g;
    $nwk =~ s/;$//;

    my $node_counter = 0;

    # Use a simple state-machine parser
    my @stack;  # stack of parent node names
    my $buf = '';

    sub finish_node {
        my ($label, $parent) = @_;
        $label =~ s/:\S+$//;   # strip branch length
        $label =~ s/^\s+|\s+$//g;
        return $label;
    }

    # Iterative Newick parser
    my $pos = 0;
    my $cur_parent = undef;

    while ($pos < length($nwk)) {
        my $ch = substr($nwk, $pos, 1);

        if ($ch eq '(') {
            push @stack, $cur_parent;
            $cur_parent = undef;
            $buf = '';
            $pos++;
        }
        elsif ($ch eq ')') {
            # End of children list; what follows is the internal node label
            $pos++;
            # read label + optional branch length
            my $label = '';
            while ($pos < length($nwk) && substr($nwk,$pos,1) !~ /[,);(]/) {
                $label .= substr($nwk,$pos,1);
                $pos++;
            }
            $label =~ s/:\S+$//;
            $label =~ s/^\s+|\s+$//g;
            $node_counter++ unless $label;
            my $iname = $label || "Node$node_counter";
            my $par = pop @stack;
            $parent_of{$iname} = $par if defined $par;
            push @{ $children{$par} }, $iname if defined $par;
            $cur_parent = $iname;
            $buf = '';
        }
        elsif ($ch eq ',') {
            # Finish current leaf/buffer
            if ($buf ne '') {
                my $label = $buf;
                $label =~ s/:\S+$//;
                $label =~ s/^\s+|\s+$//g;
                if ($label ne '') {
                    $parent_of{$label} = $cur_parent if defined $cur_parent;
                    push @{ $children{$cur_parent} }, $label if defined $cur_parent;
                }
            }
            $buf = '';
            $pos++;
        }
        else {
            $buf .= $ch;
            $pos++;
        }
    }
    # Handle remaining buffer (root label or lone leaf)
    if ($buf ne '') {
        my $label = $buf;
        $label =~ s/:\S+$//;
        $label =~ s/^\s+|\s+$//g;
        if ($label ne '') {
            $parent_of{$label} = $cur_parent if defined $cur_parent;
            push @{ $children{$cur_parent} }, $label if defined $cur_parent;
        }
    }
}

# Collect all node names
my %all_nodes;
$all_nodes{$_}++ for keys %parent_of;
$all_nodes{$_}++ for keys %children;

# Find root (no parent)
my ($root) = grep { !defined $parent_of{$_} || $parent_of{$_} eq '' } keys %all_nodes;

# Compute depth from root
my %depth;
{
    my @queue = ($root);
    $depth{$root} = 0;
    while (@queue) {
        my $n = shift @queue;
        for my $ch (@{ $children{$n} // [] }) {
            $depth{$ch} = ($depth{$n} // 0) + 1;
            push @queue, $ch;
        }
    }
}

# --------------------------------------------------------------------------
# Read IQ-TREE state file → internal node states {node => [base_at_site1, ...]}
# --------------------------------------------------------------------------
my %node_states;  # {node_name => [base1, base2, ...]}  (1-indexed sites → 0-indexed array)
{
    open my $fh, '<', $state_file or die "Cannot open $state_file: $!";
    while (<$fh>) {
        next if /^#/ || /^Node\tSite/;
        chomp;
        my @cols = split /\t/, $_;
        next unless @cols >= 3;
        my ($node, $site, $state) = @cols[0,1,2];
        next unless $site =~ /^\d+$/ && $node =~ /\S/;
        $node_states{$node}[$site - 1] = uc($state);
    }
    close $fh;
}

# Also add leaf states from the alignment
for my $leaf (keys %leaf_seq) {
    my $seq = $leaf_seq{$leaf};
    for my $i (0 .. length($seq) - 1) {
        $node_states{$leaf}[$i] = uc(substr($seq, $i, 1));
    }
}

# --------------------------------------------------------------------------
# Extract per-branch mutations
# --------------------------------------------------------------------------
# Collect all mutations: {mut_str => [node_names]}  for homoplasy detection
my %mut_nodes;

# Output lines for .db file
my @db_lines;

for my $node (sort keys %all_nodes) {
    my $par = $parent_of{$node};
    next unless defined $par && $par ne '';

    my $node_seq = $node_states{$node}   // [];
    my $par_seq  = $node_states{$par}    // [];

    my $d = $depth{$node} // 0;

    for my $i (0 .. $n_sites - 1) {
        my $nb = $node_seq->[$i] // 'N';
        my $pb = $par_seq->[$i]  // 'N';
        next if $nb eq 'N' || $pb eq 'N' || $nb eq '?' || $pb eq '?';
        next if $nb eq $pb;

        my $gpos  = $positions[$i];
        my $mut   = "${gpos}_${nb}";
        push @db_lines, "$node\t$d\t$mut\tdb";
        push @{ $mut_nodes{$mut} }, $node;
    }
}

# --------------------------------------------------------------------------
# Write .db file
# --------------------------------------------------------------------------
open my $fh_db, '>', $out_db or die "Cannot write $out_db: $!";
print $fh_db "$_\n" for @db_lines;
close $fh_db;

# --------------------------------------------------------------------------
# Write homoplasy file (mutations found on >=2 independent branches)
# --------------------------------------------------------------------------
open my $fh_homo, '>', $out_homo or die "Cannot write $out_homo: $!";
for my $mut (sort keys %mut_nodes) {
    my @nodes = @{ $mut_nodes{$mut} };
    if (@nodes >= 2) {
        print $fh_homo join("\t", $mut, scalar(@nodes), join(',', @nodes)), "\n";
    }
}
close $fh_homo;
