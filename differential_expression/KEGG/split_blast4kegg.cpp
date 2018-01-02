#include <iostream>
#include <string>
#include <fstream>
#include <map>
#include <vector>
#include <getopt.h>
#include <sys/stat.h>
#include <unistd.h>
#include <stdio.h>
#include <dirent.h>
#include <stdlib.h>
#include <boost/regex.hpp>
#include <boost/algorithm/string_regex.hpp>
using namespace std;
int Usage(){
	cerr << "Usage: split_blast4kegg -b blast_result -o output_directory" << endl;
	return 0;
}
vector<string> SplitStr(std::string str, char delim){
	std::vector<std::string> vs;
	std::string temp;
	std::string new_str = str;
	for(int i=0; i<new_str.length(); i++){
		if(new_str[i]!=delim){
			temp+=new_str[i];
		}else{
			if(temp.compare("") != 0){ vs.push_back(temp);}
			temp = "";
		}
	}
	if(temp.compare("") != 0){ vs.push_back(temp);}
	return vs;
}

int main(int argc, char *argv[]){
	optind = 0;
	const char * shortOpt = "b:o:h";
	struct option longOpt[]={
	    {"blast", 1, NULL, 'b'},
		{"odir", 1, NULL, 'o'},
		{"help", 0, NULL, 'h'},
	};
	string blast("");
	string odir("");
	int nextOpt;
	while ((nextOpt = getopt_long(argc, argv, shortOpt, longOpt, NULL)) != -1){
		switch (nextOpt){
			case 'b':
				blast = optarg;
				break;
			case 'o':
				odir  = optarg;
				break;
			case 'h':
				Usage();
		}
	}
	if(blast == "" || odir == ""){
		return Usage();
	}
	ifstream hblast(blast.c_str());
	string line;
	map<string, vector<string> > all_info;
	map<string, vector<string> >::iterator it;
	while(getline(hblast, line)){
		vector<string> rs;
		boost::split( rs, line, boost::is_any_of( "\t" ), boost::token_compress_on );
		vector<string> temp=SplitStr(rs[1], ':');
		it = all_info.find(temp[0]);
		all_info[temp[0]].push_back(line);
		//cout << temp[0] << "\t" << line << endl;
		//if(it==all_info.end()){
		//	all_info[temp[0]] = line;
		//}else{
		//	all_info[temp[0]] = all_info[temp[0]]+line;
		//}
	}
	hblast.close();
	for(map<string, vector<string> >::iterator ti=all_info.begin(); ti!=all_info.end(); ti++){
		string ofilename = odir+"/"+ti->first+"_blast.txt";
		ofstream ofile(ofilename.c_str());
		vector<string>::iterator mapvec_itor = ti->second.begin();
		for(;mapvec_itor!=ti->second.end(); mapvec_itor++){
			ofile << *mapvec_itor << endl;
		}
		ofile.close();
		//cout << "export PYTHONPATH=/share/work3/minghb/software/kobas-3.0/src:$PYTHONPATH && python " << bin << "/annotate.py -i " << odir << "/" << ti->first << "_blast.txt" << " -t blastout:tab -s " << ti->first << " -o " << odir << "/" << ti->first << ".annotate && rm " << odir << "/" << ti->first << "_blast.txt" << endl;
	}
}
