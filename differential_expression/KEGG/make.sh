#/share/public/software/gcc-5.4/bin//g++ -o split_blast4kegg -I/share/public/software/boost1.57.0/include -L/share/public/software/gzstream -L/share/public/software/boost1.57.0/lib/ -lboost_regex split_blast4kegg.cpp
source /share/public/config/ssconfig/bashrc_PB
./split_blast4kegg  -b /share/work2/liuying/project/iso-seq/BFC2013404/result/edgeR/new.RefIso.fa.KEGG.blast.xls -o edgeR.12345.dir/KEGG/A_vs_B/total/
