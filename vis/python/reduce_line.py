#! /usr/bin/env python



f = open('list', 'r')
outf=open('newlist','w')
count=0 
for line in f:
    res = count%4 
    if res==3:
        outf.write(line)
    count=count+1

f.close()
outf.close()



