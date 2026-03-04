#!/usr/bin/perl
use strict; 
use warnings;

#remove loci at which N(>threshold) strains have ambiguous base
#by using fasta file which composed by all of the strain.

die "usage:perl $0 <cfa_file> <pos_file> <output_cfa> <output_pos>\n" if @ARGV==0;

open (FA, "<$ARGV[0]");

my $hpt = $ARGV[0];
$hpt =~ s/\.fa|.fasta|.fna//;

my $thr = 0.05;

my $output_fa = $ARGV[2];
my $output_pos = $ARGV[3];
open (HPT, ">$output_fa");
open (LOC, ">$output_pos");

my %seq;
my $str;
my $len;
while(<FA>){
	chomp;
	$_ =~ s/\s+//g;
	if ($_ =~ /^>/){
		$str = $_;
	}elsif($_){
		$seq{$str} = $_; 
		$len = length($_);
	}
}
close FA;

my %loc;
my %ref;
my $n = 0;

my $sn = keys %seq;
my @loci;
$len = $len - 1; 
for my $loc (0..$len){
	my %gt=();
	$gt{N}=0;
	foreach my $k(keys %seq){
		my $nuc = substr ($seq{$k}, $loc, 1);
        	if($nuc =~ m/A|T|G|C|a|t|g|c|n|N/){
			$nuc= uc($nuc);    #capitalise "a t g c n"
       		}
        	if($nuc =~ m/N|-|\?/){     
			$gt{N}++;
		}
		elsif(!exists $gt{$nuc}){
			$gt{$nuc}=1; 
		}else{
			$gt{$nuc}++;
		}
	}
	my $frqmis = $gt{N}/$sn;
	my $gtn = keys %gt; 
    if($gt{N} == 0 && $frqmis < $thr && $gtn >= 2){    #nucletide type should be at least 2, if only 1, no useful information provided.
        push @loci, $loc;
    }elsif($gt{N} >0 && $frqmis < $thr && $gtn > 2){
		push @loci, $loc;
	}
}

open (POS, "<$ARGV[1]");
my @pos;
while(<POS>){
	chomp;
	push @pos,$_
}

foreach my $k(keys %seq){
        my $j;
        $j=$k;
        $j=~s/\.cns//;
        print HPT $j, "\n";
        foreach (@loci){
                my $nuc = substr ($seq{$k},$_,1);
                print HPT $nuc;
        }
        print HPT "\n";
}

print scalar(@loci)." check loci length\n";

foreach(@loci){
    print LOC $pos[$_],"\n";
}

close HPT;
close LOC;
