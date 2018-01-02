import os,sys
f1=open(sys.argv[1],'r')
f2=open(sys.argv[2],'r')
dict={}
result=[]
title = f2.next().rstrip() + "\tssr_position\n"
result.append(title)
for each1 in f1:
	each1=each1.strip().split("\t")
	dict[each1[0]]='\t'.join(each1[1::])
for each2 in f2:
	each2=each2.strip().split("\t")
	if each2[0] in dict.keys():
		dict1=dict[each2[0]].split('\t')
		if int(each2[5])>=min(int(dict1[0]),int(dict1[1])) and int(each2[6])<=max(int(dict1[0]),int(dict1[1])):
			text='\t'.join(each2[0::])+'\t'+'utr3'
		elif int(each2[5])>=min(int(dict1[2]),int(dict1[3])) and int(each2[6])<=max(int(dict1[2]),int(dict1[3])):
			text='\t'.join(each2[0::])+'\t'+'utr5'
		elif int(each2[5])>=min(int(dict1[4]),int(dict1[5])) and int(each2[6])<=max(int(dict1[4]),int(dict1[5])):
			text='\t'.join(each2[0::])+'\t'+'cds'
		else:
			text='\t'.join(each2[0::])+'\t'+'undetermined'
	else:
		text='\t'.join(each2[0::])+'\t'+'undetermined'
        result.append(text+'\n')
open(sys.argv[3],'w').writelines(result)
