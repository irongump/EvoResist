#!/usr/bin/perl
# db2mutation.pl
#
# Convert the .db file (branch mutation database) produced by
# nodes_base_locus_iqtree.pl to the standard mutation text format.
#
# Usage:
#   perl db2mutation.pl <input.db>
#
# The .db file has tab-separated columns:
#   node_id  depth  POS_ALT  db
#
# Output (stdout) is the same content, printed line by line.

use strict;
use warnings;

my $db_file = $ARGV[0] or die "Usage: $0 <input.db>\n";

open my $fh, '<', $db_file or die "Cannot open '$db_file': $!\n";
while (<$fh>) {
    print;
}
close $fh;
