#!usr/bin/perl
use warnings;

my $strand;

open F1, $ARGV[0] or die $!;
while(<F1>){
	chomp;
	@a=split " ",$_;
	if($#a==3){
		$a[1]=~s/\[gene=//;
		$a[1]=~s/\]//;
		$a[2]=~s/\[locus_tag=//;
		$a[2]=~s/\]//;
		$a[3]=~s/\[location=//;
		$a[3]=~s/\]//;
		if($a[3]=~m/complement/){
			$strand="-";
			$a[3]=~s/complement//;
			$a[3]=~s/\(//;
			$a[3]=~s/\)//;
#			print "$a[3]\n";
			$a[3]=~s/\.\./_/;
			@b=split "_",$a[3];
			}else{
				$strand="+";
				$a[3]=~s/\(//;
				$a[3]=~s/\)//;
				$a[3]=~s/\.\./_/;
				@b=split "_",$a[3];
				}
		print "$a[2]\t$a[1]\t$strand\t$b[0]\t$b[1]\n";				
		}elsif($#a==2){
			$a[1]=~s/\[gene=//;
			$a[1]=~s/\[locus_tag=//;
			$a[1]=~s/\]//;
			$a[2]=~s/\[location=//;
			$a[2]=~s/\]//;
			if($a[2]=~m/complement/){
			$strand="-";
			$a[2]=~s/complement//;
			$a[2]=~s/\(//;
			$a[2]=~s/\)//;
#			print "$a[3]\n";
			$a[2]=~s/\.\./_/;
			@b=split "_",$a[2];
			}else{
				$strand="+";
				$a[2]=~s/\(//;
				$a[2]=~s/\)//;
				$a[2]=~s/\.\./_/;
				@b=split "_",$a[2];
				}
		print "$a[1]\t-\t$strand\t$b[0]\t$b[1]\n";
			}
		
	}
close F1;	