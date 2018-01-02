import sys
file1=open(sys.argv[1],'r')
file2=open(sys.argv[2],'r')
file3=open(sys.argv[3],'r')
dict1={}
dict2={}
dict3={}
result=[]
for each in file1:
	each=each.strip().split('\t')
	dict1[each[0]]='\t'.join(each[1::])
for each in file2:
        each=each.strip().split('\t')
        dict2[each[0]]='\t'.join(each[1::])
for each in file3:
        each=each.strip().split('\t')
        dict3[each[0]]='\t'.join(each[1::])
for key2 in dict2:
	if key2 not in dict1.keys():
		dict1[key2]='0'+'\t'+'0'+'\t'+dict2[key2]
	else:
		dict1[key2]=dict1[key2]+'\t'+dict2[key2]
for key1 in dict1:
	if key1 not in dict2.keys():
		dict1[key1]=dict1[key1]+'\t'+'0'+'\t'+'0'
for key3 in dict3:
        if key3 not in dict1.keys():
                dict1[key3]='0'+'\t'+'0'+'\t'+'0'+'\t'+'0'+'\t'+dict3[key3]
        else:
                dict1[key3]=dict1[key3]+'\t'+dict3[key3]
for key1 in dict1:
        if key1 not in dict3.keys():
                dict1[key1]=dict3[key1]+'\t'+'0'+'\t'+'0'

for key in dict1:
        dict1[key]=key+'\t'+dict1[key]
        result+=dict1[key]+"\n"
open(sys.argv[4],'w').writelines(result)
