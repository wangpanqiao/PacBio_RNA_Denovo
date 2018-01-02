import sys
import os
import os.path
import re
import argparse
root_dir=os.getcwd()
parser = argparse.ArgumentParser(description="Prokaryotes pipline v2.0")
parser.add_argument('--project',help='project name, maybe same with the name of root dir, [REQUIRED]',required=True)
parser.add_argument('--root_dir',help='project dir',default=root_dir)
argv = vars(parser.parse_args())
project=argv['project'].strip()
root_dir=argv['root_dir'].strip()
os.chdir(root_dir+'/SSR/')

f = open('SSR_utr.sh','w')

f.write('perl /PUBLIC/source/RNA/noRef/SSR_Primer/bin/get_UTR.pl -info %s -seq %s -out %s' %(root_dir+'/ANNOTATION_ALL/CDSprediction/'+project+'.orf.info',root_dir+'/TRINITY/assembly_INFO/unigene.fasta',root_dir+'/SSR/SSR_Primer/tmp\n'))
f.write('python /PUBLIC/source/RNA/noRef/SSR_Primer/bin/exact_pos.py %s %s' %(root_dir+'/SSR/SSR_Primer/tmp/'+'unigene.3_UTR.fasta',root_dir+'/SSR/SSR_Primer/tmp/'+'3_UTR\n'))
f.write('python /PUBLIC/source/RNA/noRef/SSR_Primer/bin/exact_pos.py %s %s' %(root_dir+'/SSR/SSR_Primer/tmp/'+'unigene.5_UTR.fasta',root_dir+'/SSR/SSR_Primer/tmp/'+'5_UTR\n'))
f.write('python /PUBLIC/source/RNA/noRef/SSR_Primer/bin/exact_pos.py %s %s' %(root_dir+'/SSR/SSR_Primer/tmp/'+'unigene.cds.fasta',root_dir+'/SSR/SSR_Primer/tmp/'+'cds\n'))
f.write('python /PUBLIC/source/RNA/noRef/SSR_Primer/bin/multi_length.py %s %s %s %s' %(root_dir+'/SSR/SSR_Primer/tmp/'+'3_UTR',root_dir+'/SSR/SSR_Primer/tmp/'+'5_UTR',root_dir+'/SSR/SSR_Primer/tmp/'+'cds',root_dir+'/SSR/SSR_Primer/tmp/'+'combine.xls\n'))
f.write('python /PUBLIC/source/RNA/noRef/SSR_Primer/bin/ssr_pos.py %s %s %s' %(root_dir+'/SSR/SSR_Primer/tmp/'+'combine.xls',root_dir+'/SSR/SSR_Primer/unigene.fasta.misa',root_dir+'/SSR/SSR_Primer/'+project+'.SSR.pos.misa\n'))
f.write('rm -r %s'%(root_dir+'/SSR/SSR_Primer/tmp/'+'\n'))
f.close()

