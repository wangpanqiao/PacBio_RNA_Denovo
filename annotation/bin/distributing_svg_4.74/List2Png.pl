#!/usr/local/bin/perl -w
# 
# Copyright (c) BMK 2009
# Writer:         Yangsh <yangsh@biomarker.com.cn>
# Program Date:   2010.
# Modifier:       Yangsh <yangsh@biomarker.com.cn>
# Last Modified:  2010.
my $ver="1.0";
my $BEGIN=time();

use strict;
use Cwd;
use Getopt::Long;
use Data::Dumper;
use FindBin qw($Bin $Script);
use File::Basename qw(basename dirname);

######################请在写程序之前，一定写明时间、程序用途、参数说明；每次修改程序时，也请做好注释工作

if (@ARGV!=1) 
{
	print "Usage: list \n";
	exit;
}

###############Time
my $Time_Start;
$Time_Start = sub_format_datetime(localtime(time()));
print "\nStart Time :[$Time_Start]\n\n";
################
my $programe_dir=basename($0);
my $path=dirname($0);
my $dis_svg="$path/distributing_svg.pl";
my $svg2png="$path/svg2xxx_release/svg2xxx";

`perl $dis_svg $ARGV[0] $ARGV[0].svg `;
`perl $svg2png $ARGV[0].svg`;
###############Time
my $Time_End;
$Time_End = sub_format_datetime(localtime(time()));
print "\nEnd Time :[$Time_End]\n\n";
&Runtime($BEGIN);

#+---------------------
#        Subs         |
#+---------------------

sub Runtime
{ # &Runtime($BEGIN);
	my ($t1)=@_;
	my $t=time()-$t1;
	print "Total elapsed time: ${t}s\n";
}

sub sub_format_datetime {#Time calculation subroutine
    my($sec, $min, $hour, $day, $mon, $year, $wday, $yday, $isdst) = @_;
	$wday = $yday = $isdst = 0;
    sprintf("%4d-%02d-%02d %02d:%02d:%02d", $year+1900, $mon+1, $day, $hour, $min, $sec);
}
