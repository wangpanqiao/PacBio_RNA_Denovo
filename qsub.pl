#!/usr/bin/perl -w
use strict;
use Getopt::Long;
use File::Basename;

my $usage =
"
Usage:

  Options:
    -inputdir      <dir>         The path of project directory
    -sample        <str>         Sample name
    -fileno        <num>         The number of required input files
    -help                        Help message
    -step          <num>         which step to run from,default is 1 [option: 1,2,3,4,5,6,7,8,9,10,11,11]
                                    1: IsoSeq,
		                    2: LSC,
				    3: CD-hit-est,
                                    4: Annotation,
				    5: ORF transdecoder,
				    6: OrthoMCL,
				    7: SSR,
				    8: LncRNA,
				    9: Get reference isoform,
				   10: Expression evaluation,
				   11: Differential expression,
                                   12: Package result and web report
    -cut           <num>         which step to cut-off end running,default is 12 [option: 1,2,3,4,5,6,7,8,9,10,11,12]
For example:
    perl $0 -inputdir /home/xingsk/workdir/project/PacBio_RNA_Denovo/project_test -sample BR20160620 -step 1 -cut 12
";
my ($shell_dir,$sample,$fileno,$step,$cut,$help,$finish,$cmdout,$sh,$mark);
GetOptions(
                                  "inputdir=s"=>\$shell_dir,
                                  "sample=s"=>\$sample,
                                  "fileno=i"=>\$fileno,
                                  "step=i"=>\$step,
                                  "cut=i"=>\$cut,
                                  "h|help"=>\$help
);
$step ||= 1;
$cut ||= 12;
if (!$shell_dir || !$sample || $help){
	die "$usage\n";
}
die "start-step must bigger than the end-step\n" if ($step > $cut);

#####################################################################################
#                     check current running status                                  #
#####################################################################################

if(-e "$shell_dir/worksh/failed")
{exit(0);}
#####################################################################################
#                            qsub IsoSeq.sh                                         #
#####################################################################################
my @names=split /\//,$sample;
if ($step<2 && $cut>0){
         system("cd $shell_dir/result/IsoSeq/$sample && qsub -cwd -l vf=15G IsoSeq.sh");
}
#####################################################################################
#                              qsub LSC.sh                                          #
#####################################################################################
if ($step<3 && $cut>1){
        $finish = 1;
	while($finish){
		sleep 60;
		$cmdout = `ls -l $shell_dir/worksh/ | grep -c '$sample.IsoSeq.mark\$'` or die "$!\n";
		if ($cmdout == $fileno){
			$finish = 0;
		}
	}
	system("cd $shell_dir/result/LSC/$sample && qsub -cwd LSC.sh");
}

#####################################################################################
#                            qsub CD-hit-est.sh                                     #
#####################################################################################

if ($step<4 && $cut>2){
	unless (-e "$shell_dir/result/LSC/$sample"){
        	$finish = 1;
		while($finish){
		sleep 60;
		$cmdout = `ls -l $shell_dir/worksh/ | grep -c '$sample.IsoSeq.mark\$'` or die "$!\n";
		if ($cmdout == $fileno){
			$finish = 0;
					}
             			}
		}	
	else {
        	$finish = 1;
		while($finish){
		sleep 60;
		if (-e "$shell_dir/result/LSC/$sample/LSC.mark"){ $finish = 0;}
  		}
	     }	
	system("cd $shell_dir/result/CD-hit-est/$sample  && qsub -cwd CD-hit-est.sh");
}


#####################################################################################
#                             qsub Annotation.sh                                    #
#####################################################################################
if ($step<5 && $cut>3){
	$finish = 1;
	while($finish){
		sleep 60;
		if (-e "$shell_dir/result/CD-hit-est/$sample/CD-hit-est.mark"){
			$finish = 0;
		}
	}
	system("cd $shell_dir/result/Annotation/$sample && qsub -cwd -l rf=15G Annotation.sh");
}

#####################################################################################
#                            qsub Transdecoder.sh                                   #
#####################################################################################
if ($step<6 && $cut>4){
	$finish = 1;
	while($finish){
		sleep 60;
		if (-e "$shell_dir/result/CD-hit-est/$sample/CD-hit-est.mark"){
			$finish = 0;
		}
	}
		system("cd $shell_dir/result/Transdecoder/$sample && qsub -cwd Transdecoder.sh");
}

#####################################################################################
#                              qsub OrthoMCL.sh                                     #
#####################################################################################
if ($step<7 && $cut>5){
	$finish = 1;
if($fileno==1){
while($finish){
                sleep 60;
                if (-e "$shell_dir/result/Transdecoder/$sample/Transdecoder.mark"){
                        $finish = 0;
                }
        }
                system("cd $shell_dir/result/OrthoMCL/$sample && sh OrthoMCL.sh");##OrthoMCL.$names[1].sh contain "make" and "qsub" actions
}else{
        $fileno=$fileno-1;
while($finish){
		sleep 60;
		if (-e "$shell_dir/worksh/OrthoMCL.mark.$fileno"){
			$finish = 0;
		}
	}
		system("cd $shell_dir/result/OrthoMCL/$sample && sh OrthoMCL.sh");##OrthoMCL.sh contain "make" and "qsub" actions
}
}
#####################################################################################
#                                 qsub SSR.sh                                       #
#####################################################################################
if ($step<8 && $cut>6){
        $finish = 1;
        while($finish){
                sleep 60;
                if (-e "$shell_dir/result/CD-hit-est/$sample/CD-hit-est.mark"){
                        $finish = 0;
                }
        }
                system("cd $shell_dir/result/SSR/$sample && qsub -cwd SSR.sh");
}

#####################################################################################
#                                qsub LncRNA.sh                                     #
#####################################################################################
if ($step<9 && $cut>7){
        $finish = 1;
        while($finish){
                sleep 60;
                if (-e "$shell_dir/result/CD-hit-est/$sample/CD-hit-est.mark"){
                        $finish = 0;
                }
        }
                system("cd $shell_dir/result/LncRNA/$sample && qsub -cwd LncRNA.sh");
}
#####################################################################################
#                               qsub RefIso.sh                                      #
#####################################################################################
if ($step<10 && $cut>8){
        $finish = 1;
        while($finish){
                sleep 60;
                $cmdout = `ls -l $shell_dir/worksh/ | grep -c '.CD-hit-est.mark\$'` or die "$!\n";
                if ($cmdout == $fileno){
                        $finish = 0;
                }
  }
        system("cd $shell_dir/result/RefIso/  && qsub -cwd RefIso.sh");
}
#####################################################################################
#            qsub rsem-prepare-reference.sh and Sample_expression.sh                #
#####################################################################################
if ($step<11 && $cut>9){
	$finish = 1;
unless (-e "$shell_dir/result/RSEM/reference_prepare_start.mark"){
       	while($finish){
		if (-e "$shell_dir/result/RefIso/RefIso.mark" ){
			$finish = 0;}
                sleep 1;
	}
		system("cd $shell_dir/result/RSEM  && qsub -cwd rsem-prepare-reference.sh");
	}
 
	$finish = 1;
       	while($finish){
		sleep 60;
		if (-e "$shell_dir/result/RSEM/reference_prepare.mark"){
			$finish = 0;
		}
	}
		system("cd $shell_dir/result/RSEM  && qsub -cwd expression.$sample.sh");
	}
 
#####################################################################################
#                               qsub diff_exp.sh                                    #
#####################################################################################
if ($step<12 && $cut>10){
	$finish = 1;
	while($finish){
		sleep 60;
		$cmdout = `ls -l $shell_dir/result/RSEM/ | grep -c 'expression.mark\$'` or die "$!\n";
		if ($cmdout == $fileno){
			$finish = 0;
		}
		
	}
		system("cd $shell_dir/result/edgeR  && qsub -cwd diff_exp.sh");
}

#####################################################################################
#                              qsub Web_report.sh                                   #
#####################################################################################
if ($step<13 && $cut>11){
	$finish = 1;
	while($finish){
		sleep 60;
		$sh = `ls -l $shell_dir/worksh/ | grep -c '.sh\$'` or die "$!\n";
		$mark = `ls -l $shell_dir/worksh/ | grep -c '.mark\$'` or die "$!\n";
		if($sh == $mark){	
                       $finish = 0;
		}
	}
		system("cd $shell_dir  && qsub -cwd -l h=compute-0-0 Web_report.sh");
}

 
