#! /urs/bin/perl -w
use strict;
use File::Basename qw(basename dirname);

my (%hash,%pathway);
my @file =glob "$ARGV[0]/*.annotate";
foreach my $file (@file){
	open FILE,"$file";
	$/="--------------------";<FILE>;$/="////";<FILE>;
	while (<FILE>){
		my @line =split /\n/,$_;
		if (@line >5){
			my ($gene,$pathway,$info);
			foreach my $item(@line){
				$item=~s/\n//g;
				my @con=split /\t/,$item;
				if (defined $con[0]){
					if ($con[0]=~/^Query/){
						$gene=$con[1];
					}
					if ($con[0]=~/^Pathway/ ||$con[0]=~/^\s+/){
						$pathway=$con[1];
						$info="$gene\t$pathway";
						$hash{$info}=1;
					}
				}
			}
		}
	}
	close FILE; 
}
open ALL,"$ARGV[1]";
$/="--------------------";
my $head =<ALL>;
print $head;
$/="////";
my $head2 = <ALL>;
print "$head2\n";

while (<ALL>){
	my @line =split /\n/,$_;
	if (@line < 3){
		print "$_";
	}else{  
		my ($gene,$pathway,$info);
		foreach my $item(@line){
			$item=~s/\n//;
			my @con=split /\t/,$item;
			if (defined $con[0]){
				if ($con[0]=~/^Query/){
					$gene=$con[1];
					print "$item\n";
				}elsif($con[0]=~/^Pathway/ || $con[0] =~ /^\s+/){
					
					$pathway=$con[1];
					$info="$gene\t$pathway";
					if (exists $hash{$info}){
						if (exists $pathway{$gene}){
							shift @con;
							my $cont=join "\t",@con;
							print "                    \t$cont\n";		
						}else{
							shift  @con;
							my $cont=join "\t",@con;
							print "Pathway:            \t$cont\n";
							$pathway{$gene}=1;
						}
					}
				}else {print "$item\n";}
			} 
		}
	}	
}
close ALL;
