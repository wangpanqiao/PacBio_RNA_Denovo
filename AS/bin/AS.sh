python /share/work1/staff/chaijc/software/IsoSeq_AS_de_novo/1step.py all_vs_all_blast.xml 1.xls
python /share/work1/staff/chaijc/software/IsoSeq_AS_de_novo/2.1step.py 1.xls 2.1.xls
python /share/work1/staff/chaijc/software/IsoSeq_AS_de_novo/2.2step.py 2.1.xls 2.2.xls
python /share/work1/staff/chaijc/software/IsoSeq_AS_de_novo/3step.py 2.2.xls 3.xls
python /share/work1/staff/chaijc/software/IsoSeq_AS_de_novo/4step.py 3.xls 4.xls
python /share/work1/staff/chaijc/software/IsoSeq_AS_de_novo/5step.py 4.xls 5.xls
python /share/work1/staff/chaijc/software/IsoSeq_AS_de_novo/6.1step.py 5.xls 6.1.xls
sed "s/No definition line//g" 6.1.xls > 6.1step.clean.txt
sort -k 2,2V -k 3,3rn 6.1step.clean.txt > 6.1step.clean.sortS.txt
awk '{if(! a[$2]){print; a[$2]++}}' 6.1step.clean.sortS.txt > 6.1step.clean.sortS.uniqS.txt
sort -k 1,1V -k 3,3rn 6.1step.clean.sortS.uniqS.txt > 6.1step.clean.sortS.uniqS.sortQ.txt
awk '{if(! a[$1]){print; a[$1]++}}' 6.1step.clean.sortS.uniqS.sortQ.txt > 6.1step.clean.sortS.uniqS.sortQ.uniqQ.txt
