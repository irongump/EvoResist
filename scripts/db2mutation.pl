#Filter host selection phylogeny nodes and SNVs.
#Average SNVs counts from node to tips were calculated,
#If the counts > e.g. 0.5 SNV/year * (2021-1940) = 40:
#   Output the node and SNVs
#else:
#   ignore the node

use warnings;
use strict;
use List::Util qw/max min/;

die "usage:perl $0\n" if @ARGV==0;

open DB,"<$ARGV[0]" or die "db file $!\n";
$/ = '>';
my %node_count; #node to snv counts
my %node_snvs; #node to snvs
<DB>;
while (<DB>){
    chomp;
    my @lines = split /\n/;
    my $nodeid = $lines[0];
    if ($lines[1] =~/no site/){
        $node_count{$nodeid} = 0;
        $node_snvs{$nodeid} = 0;
    }
    else {
        my @sites = split /\s+/,$lines[1];
        my @snvs = split //,$lines[5];
        my $num = scalar @sites;
        my @mut; #mutations
            for (0..$#sites){
                push @mut, $sites[$_]."_".$snvs[$_];
            }
        $node_count{$nodeid} = $num;
        $node_snvs{$nodeid} = \@mut;
    }
}
close DB;

for (keys %node_snvs){
    my $node = $_;
    if ($node_count{$node} !=0){ 
        for (@{$node_snvs{$node}}){
            print "$node\t$node_count{$node}\t$_\tdb\n";
        }
    }
}
