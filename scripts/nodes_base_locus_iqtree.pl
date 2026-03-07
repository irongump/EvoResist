#!~/.conda/envs/bioperl/bin/perl
use warnings;
use strict;
use Bio::TreeIO;
use List::Util qw/max min/;

die "usage:perl $0 <iqtree_tree_file> <locus_file> 
<ancestor_seq.file> <raw.fasta> <output_db_file> <output_homoplasy_file>\n" if @ARGV==0;

my $treefile = $ARGV[0];
my $treeio = new Bio::TreeIO(-format=>'newick',-file=>$treefile);
my (%ancestor_descendant, %strain_node, %node_sequence, %hash_all, %hash_descendant, $total_node);
if (my $tree = $treeio->next_tree) {
    my @nodes = $tree->get_nodes;
    for (@nodes){
        my $des = $_;
        my $anc = $des->ancestor;
        my $outdes = $des->id;
        if (defined $anc){
            my $outanc = $anc->id;
            if ($outanc =~ /.+/ and $outdes =~ /.+/){
                print "$outanc\t$outdes\n";
                push @{$ancestor_descendant{$outanc}}, $outdes;
                $hash_descendant{$outdes} = $outanc;
                $hash_all{$outanc} = 1;
                $hash_all{$outdes} = 2;
            }
        }
    }
    $total_node = @nodes;
}

#my @desset = ('Node1');
#$ancestor_descendant{'Node0'} = \@desset;

my (@node_numeric, @node_text);
for (keys %hash_all){
    my $node_type = $_;
    if ($node_type =~ /^Node\d+$/i) {
        $node_type =~ s/Node//i;
        push @node_numeric, $node_type;
    }
    else {
        push @node_text, $node_type;    
    }
}

my $max_node = (max @node_numeric);
#my %number_node;
#map {$number_node{'node'.$_} = 1} @node_numeric;

my $node_text_count = @node_text;
#for (0..$#node_text){
#    my $i = $_;
#    my $strain = $node_text[$i];
#    for (1..$total_node){
#        my $j = $_;
#        if (!exists $number_node{'node'.$j}){
#            $strain_node{$strain} = 'node'.$j;
#            $number_node{'node'.$j} = 1;
#            last;
#        }      
#    }
#}
#my $strain_node_count = scalar keys %strain_node;
#print "max:$max_node\t$node_text_count\t$strain_node_count\n";
#die "strain_node_count ne node_text_count!\n" if $node_text_count != $strain_node_count; 

my @node_n;
open ANS,"<$ARGV[2]" or die "ancestor seq file $!\n";
while (<ANS>){
    chomp;
    if (/Node/ and $_ !~ /Site/){
        my @node = split;
        #push @node_n, $node[0];
        $node_sequence{$node[0]} .= $node[2]
    
    }
}
close ANS;
my $node_count_N = scalar keys %node_sequence;#@node_n;
my $node_seq_count = scalar keys %node_sequence;
#die "sequence count:$node_seq_count ne node count:$node_count_N!\n" if $node_count_N != $node_seq_count; 

open FA,"<$ARGV[3]" or die "fa $!\n"; # fasta file of last phylogeny construct step
while (<FA>){
    chomp;
    my $name = $_;
    $name =~ s/^>//;
    #print "fa:$name\n";
    my $seq = <FA>;
    chomp $seq;
    #if ($name =~ /tb.ancestor/){
    #    $node_sequence{"Node0"} = $seq;
    #}
    #else{
    $node_sequence{$name} = $seq;
    #}
}
close FA;


open LOC,"<$ARGV[1]" or die "$!\n";
my @locus;
while(<LOC>){
    chomp;
    my @array = split;
    push @locus, $array[0];
}
close LOC;

for (keys %hash_descendant){
    my $descendant = $_;
    if (!exists $node_sequence{$descendant}){
        $node_sequence{$descendant} = $node_sequence{$hash_descendant{$descendant}};    
    }
}


my (%node_specifical_sequence, %base_with_multiple_node,%node_include_base_with_multiple_node);
for (keys %ancestor_descendant){
    my $ancestor_node = $_;
    my $ancestor_seq;
    if (exists $node_sequence{$ancestor_node}){
        $ancestor_seq = $node_sequence{$ancestor_node};
    } else{
        die "$ancestor_node:ancestor node not found!\n";
    }
    for (@{$ancestor_descendant{$ancestor_node}}){
        my $descendant_node = $_;
        my $descendant_seq = (exists $node_sequence{$descendant_node}) ? $node_sequence{$descendant_node} : $ancestor_seq;
        my $seq_length = length $descendant_seq;
        my $anc_seq_length = length $ancestor_seq;
        print "$descendant_node--\t$anc_seq_length\t$seq_length\t$#locus\n";
        die "$descendant_node:site number do not equal sequence length!\n" if ($seq_length != $#locus + 1);
        my $array_number = $seq_length - 1;
        for (0..$array_number){
            my $i = $_;
            my $ancestor_base = substr($ancestor_seq, $i, 1);
            my $descendant_base = substr($descendant_seq, $i, 1);
            if (($ancestor_base ne $descendant_base) and ($descendant_base !~
                /[N?-]/)){
                my $site = $locus[$i];
                push @{$base_with_multiple_node{"$site"."_$descendant_base"}}, $descendant_node;    
                ${$node_specifical_sequence{$descendant_node}}{$site} = $descendant_base;    
            }
        }
    }
}

my $base_with_multiple_node_count = 0;
open HOMO, ">$ARGV[5]" or die "$!\n";
print HOMO "homoplasy bases:\n";
for (keys %base_with_multiple_node){
    my $key = $_;
    my $j = @{$base_with_multiple_node{$key}};
    if ($j>1){
        $node_include_base_with_multiple_node{$key} = 1;
        $base_with_multiple_node_count += $j - 1;
        print HOMO "$key\t$j\t";
        map {print HOMO "$_\t"} @{$base_with_multiple_node{$key}};
        print HOMO "\n";
    }
}
close HOMO;

my $base_with_multiple_node_count_N = scalar keys %node_include_base_with_multiple_node;
print "base with multiple node:$base_with_multiple_node_count_N\n";
open DB,">$ARGV[4]" or die "$!\n";

for (keys %hash_all){
    my $pre_node = $_;  
    #my $aft_node = (exists $strain_node{$pre_node}) ? $strain_node{$pre_node} : $pre_node;
    my @array_pre_node = ();
    @array_pre_node = &Minlevel ($pre_node, %hash_descendant);
    if (@array_pre_node){
        my @reverse_array_pre_node =reverse @array_pre_node;
        pop @reverse_array_pre_node;
        my $rapn_count = @reverse_array_pre_node;
        print DB ">$pre_node\n";
        if (exists $node_specifical_sequence{$pre_node}){
            for (sort {$a<=>$b} keys %{$node_specifical_sequence{$pre_node}}){
                print DB "$_ ";
            }    
            print DB "\n";
        }
        else {
            print DB "no site\n";    
        }
        if (exists $hash_descendant{$pre_node}){
            print DB "$hash_descendant{$pre_node}\n";
        }
        else {
            print DB "wired\n";    
        }
        map {print DB "$_-"} @reverse_array_pre_node;
        print DB "\n";
        print DB "$rapn_count\n";
        if (exists $node_specifical_sequence{$pre_node}){
            for (sort {$a<=>$b} keys %{$node_specifical_sequence{$pre_node}}){
                print DB "${$node_specifical_sequence{$pre_node}}{$_}";
            }
            print DB "\n";
        }
        else {
            print DB "no base type\n";    
        }
    }
    else {
        print DB ">$pre_node\n";
        print DB "no site\n\n";
        print DB "root\n";
        print DB "-\n";
        print DB "0\n";
        print DB "no base type\n";
    }
}
my @node_seq;
sub Minlevel {
    my ($temp_node,%hash_temp) = @_;
    push @node_seq, $temp_node;
    if (exists $hash_temp{$temp_node}){
        my $string = $hash_temp{$temp_node};
        delete $hash_temp{$temp_node};
        &Minlevel ($string, %hash_temp);
    }
    else {
        my @node_seq_out = @node_seq;
        undef @node_seq;
        return @node_seq_out;
    }
}

