import os,sys
import re
f1=open(sys.argv[1],'r')
list = []
for line in f1:
	if line.strip().startswith('>'):
		a = re.search(r'>(\w+).*',line.strip())
		b = re.search(r'.*start:(\d+).*',line.strip())
		c = re.search(r'.*end:(\d+).*',line.strip())
		list.append([a.group(1),b.group(1),c.group(1)])

out = []
for each in list:
	text = str(each[0]) + '\t' + str(each[1]) + '\t' + str(each[2]) + '\n'
	out.append(text)
open(sys.argv[2],'w').writelines(out)
