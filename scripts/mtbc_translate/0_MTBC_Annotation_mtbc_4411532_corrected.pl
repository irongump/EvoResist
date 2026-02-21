#!usr/bin/perl
use warnings;

# this script is writted to corrected those false annotation for those linked mutations;

$name="$ARGV[0]";
$name.=".temp";
# create a temperory file $.temp (will be deleted in the end)
`perl /proj/qliulab/MTB_phy_db/src/mtbc_translate/1_MTBC_Annotation_mtbc_4411532.pl $ARGV[0] > $name`;
# read in the coding system
my %code;
open F0, "/proj/qliulab/MTB_phy_db/src/mtbc_translate/3_genetic_codes" or die $!;
while(<F0>){
        chomp;
        @a=split "\t",$_;
        $code{$a[0]}=$a[1];
        }
close F0;

# read in the initial translated file
my %tr;
my @loca;
open F1, "$name" or die $!;
while(<F1>){
chomp;
@a=split "\t",$_;
push @loca, $a[0];
$tr{$a[0]}=$_;
}
close F1;

# start to correct and output the corrected version
my %exis;
my %hash;

foreach $i(0..$#loca){
$j=$loca[$i];
@b=split "\t",$tr{$j};
$code=$b[7].$b[3];
if(!exists $hash{$code} && $b[3] ne "-"){
$hash{$code}=$b[0];
$exis{$b[0]}=0;
}elsif(exists $hash{$code} && $b[3] ne "-"){
$lc1=$hash{$code};
$lc2=$b[0];
$exis{$lc1}=1;
$exis{$b[0]}=1;
$hash{$code}.="-".$b[0];
$new=$hash{$code};
push @loca, $new;
$exis{$new}=0;
$seq1=$tr{$lc1};
$seq2=$tr{$lc2};
#print "$seq1\n$seq2\n";
@b1=split "\t",$seq1;
@b2=split "\t",$seq2;
@c1=split "-", $b1[5];
@c2=split "-", $b2[5];
foreach $k(0..2){
$wt=substr $c1[0],$k,1;
$mt1=substr $c1[1],$k,1;
$mt2=substr $c2[1],$k,1;
if($wt eq $mt1 && $wt eq $mt2){
$coden.=$wt;
}elsif($wt ne $mt1 && $wt eq $mt2){
$coden.=$mt1;
$wt3.=$wt;
$mt3.=$mt1;
}elsif($wt eq $mt1 && $wt ne $mt2){
$coden.=$mt2;
$wt3.=$wt;
$mt3.=$mt2;
}
}
$ori=$code{$c1[0]};
$mut=$code{$coden};
if($ori eq $mut){
$m="Synonymous";
}else{
$m="Nonsynonymous";
}
$seq3="$new\t$wt3\t$mt3\t$b2[3]\t$m-$ori-$mut\t$c1[0]-$coden\t$b2[6]\t$b2[7]\t$b2[8]\t$b2[9]";
$tr{$new}=$seq3;
$wt3="";
$mt3="";
$coden="";
}elsif($b[3] eq "-"){
$exis{$b[0]}=0;
}
}

foreach $n(@loca){
if($exis{$n}==0){
print "$tr{$n}\n";
}
}

`rm $name`;











=pod
open F2, $ARGV[1] or die $!;
while(<F2>){
chomp;
}
close F2;
=cut

=pod
open F2, $ARGV[2] or die $!;
while(<F3>){
chomp;
}
close F3;
=cut
