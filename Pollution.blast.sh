sed -n "1,$(cat -n /share/seq_dir1/Item1/pacbio/RD/ISO_seq_mouse/Sample_10_3-6K/D02_1/Analysis_Results/m160704_145251_42221_c100952062550000001823213106251664_s1_X0.1.subreads.fasta | grep ">" | sed -n "20001p" | awk '{print $1}')p" /share/seq_dir1/Item1/pacbio/RD/ISO_seq_mouse/Sample_10_3-6K/D02_1/Analysis_Results/m160704_145251_42221_c100952062550000001823213106251664_s1_X0.1.subreads.fasta | sed '$d'>20000line.fasta &&
/home/leiyoubing/work/software/ncbi-blast-2.2.29+/bin/blastn -evalue 1e-5 -max_target_seqs 10 -num_threads 8 -outfmt '7 qseqid sseqid evalue pident ppos length mismatch gapopen qstart qend sstart send bitscore staxids sscinames sskingdoms stitle' -db /share/work2/staff/zhouy/Database/NT/nt -query 20000line.fasta -out DM8_L2.txt &&
perl /home/jiangdezhi/workdir/script/BioInfo_Tools/blast/fast_blast/blastall_stat.pl DM8_L2.txt -o DM8_L2 && echo This-Work-is-Completed!
