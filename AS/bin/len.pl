#/usr/bin/perl -w
use strict;

my ($seq, $id) = @ARGV;
my (%hash, $id_new);
open IN, "$seq" or die "cannot open:$!";
while(<IN>){
	chomp;
	if(/^>(.*)/){
		$id_new = $1;
	}else{
		$hash{$id_new} .= $_;
	}
}
foreach (sort keys%hash){
	if($id eq $_){
		print "the length of $_ is ".length($hash{$_})."\n";
		#print "$hash{$_}";
	}
}
close IN;
	
