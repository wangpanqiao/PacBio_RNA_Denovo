#/usr/bim/perl -w
use warnings;

my($in) = @ARGV;
my $i = 1;
open IN, "$in" or die "cannot oepn:$!";
while(<IN>){
	chomp;
	if(/^>/){
		print ">transcript_$i\n";
		$i +=1;
	}else{
		print "$_\n";
	}
}
close IN;
	
	
