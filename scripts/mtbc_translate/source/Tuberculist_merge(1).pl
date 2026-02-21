#!usr/bin/perl
use warnings;

my %Rv;
my %MTB;

open F1,$ARGV[0] or die $!;
while(<F1>){
	chomp;
	@a=split "\t",$_;
	if($a[0]=~m/Rv/){
		$Rv{$a[0]}=$_;
		}elsif($a[0]=~m/MTB/){
			$MTB{$a[1]}=$_;
			}
	}
close F1;

open F2,$ARGV[1] or die $!;
while(<F2>){
	chomp;
	@b=split "\t",$_;
	if($b[0]=~m/Rv/){
		@c=split "\t",$Rv{$b[0]};
		print "$c[0]\t$c[1]\t$b[2]\t$b[3]\t$b[4]\t$c[2]\t$c[3]\n";
		}else{
			@c=split "\t",$MTB{$b[0]};
			print "$c[0]\t$c[1]\t$b[2]\t$b[3]\t$b[4]\t$c[2]\t$c[3]\n";
			}
	}
close F2;	