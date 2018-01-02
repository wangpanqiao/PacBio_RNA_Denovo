# !/usr/bin/env pyhton
# coding=utf-8
# Writer:         chaijingchao
# Program Date:   2016.08.11

import argparse,json,os
parser = argparse.ArgumentParser(description="stat the info of iso_seq")
parser.add_argument('--indir',help='the direction of project', required=True)
parser.add_argument('--sample',help='the name of sample', required=True)
parser.add_argument('--cDNA',help='the cDNA size of sample', required=True)
parser.add_argument('--outdir',help='the direction of output',required=True)
argv = vars(parser.parse_args())
i = argv["sample"]
j = argv["cDNA"]
indir = argv["indir"]
outdir = argv["outdir"]
dic = {}
cDNA = j.split(",")
out_roi = open((outdir + "/" + i + "_roi_stat.xls"),"w")
out_roi.write("Sample\tcDNA size\tReads of Insert\tRead Bases of Insert\tMean Read Length of Insert\tMean Read Quality of Insert\tMean Number of Passes\n")
out_classify = open((outdir + "/" + i + "_classify_stat.xls"),"w")
out_classify.write("Sample\tcDNA size\tNumber of reads of insert\tNumber of five prime reads\tNumber of three prime reads\tNumber of poly-A reads\tNumber of filtered short reads\tNumber of non-full-length reads\tNumber of full-length reads\tNumber of full-length non-chimeric reads\tFull-Length Percentage (FL%)\tAverage full-length non-chimeric read length\n")
out_cluster = open((outdir + "/" + i + "_cluster_stat.xls"),"w")
out_cluster.write("Sample\tcDNA size\tNumber of consensus isoforms\tAverage consensus isoforms read length\tNumber of polished high-quality isoforms\tNumber of polished low-quality isoforms\n")
for each in cDNA:
    iso_dir = indir + i + "/" + i + "_" + each + "/results"
    f1 = open(iso_dir + "/reads_of_insert_report.json")
    s1 = json.load(f1)
    roi_number_reads = s1["attributes"][0]["value"]
    roi_number_bases = s1["attributes"][1]["value"]
    roi_number_mean = s1["attributes"][2]["value"]
    roi_number_accuracy = s1["attributes"][3]["value"]
    roi_number_passes = s1["attributes"][4]["value"]
    out_roi.write(str(i) + "\t" +str(each) + "\t" + str(roi_number_reads) + "\t" + str(roi_number_bases) + "\t" + str(roi_number_mean) + "\t" + str(roi_number_accuracy) + "\t" + str(roi_number_passes) + "\n")
    f2 = open(iso_dir + "/isoseq_classify.json")
    s2 = json.load(f2)
    classify_number_reads = s2["attributes"][0]["value"]
    classify_number_five_reads = s2["attributes"][1]["value"]
    classify_number_three_reads = s2["attributes"][2]["value"]
    classify_number_poly_a_reads = s2["attributes"][3]["value"]
    classify_number_short_reads = s2["attributes"][4]["value"]
    classify_number_nfl_reads = s2["attributes"][5]["value"]
    classify_number_fl_reads = s2["attributes"][6]["value"]
    classify_number_fln_reads = s2["attributes"][7]["value"]
    classify_fln_percentage = "%.0f%%" %(float(classify_number_fln_reads)/classify_number_reads*100)
    classify_number_nfln_read_length = s2["attributes"][8]["value"]
    out_classify.write(str(i) + "\t" +str(each) + "\t" + str(classify_number_reads) + "\t" + str(classify_number_five_reads) + "\t" + str(classify_number_three_reads) + "\t" + str(classify_number_poly_a_reads) + "\t" + str(classify_number_short_reads) + "\t" + str(classify_number_nfl_reads) + "\t" + str(classify_number_fl_reads) + "\t" + str(classify_number_fln_reads) + "\t" + str(classify_fln_percentage) + "\t" + str(classify_number_nfln_read_length) + "\n")
    f3 = open(iso_dir + "/isoseq_cluster.json")
    s3 = json.load(f3)
    cluster_number_consensus = s3["attributes"][0]["value"]
    cluster_consensus_length = s3["attributes"][1]["value"]
    cluster_number_hqv = s3["attributes"][2]["value"]
    cluster_number_lqv = s3["attributes"][3]["value"]
    out_cluster.write(str(i) + "\t" +str(each) + "\t" + str(cluster_number_consensus) + "\t" + str(cluster_consensus_length) + "\t" + str(cluster_number_hqv) + "\t" + str(cluster_number_lqv) + "\n")
out_roi.close()
out_classify.close()
out_cluster.close()
