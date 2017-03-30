import os
import sys
from subprocess import *

# sf_support_index=sys.argv[1]
# sf_all_pos=sys.argv[2]
# sf_align=sys.argv[3]
#
# m_all_pos=[]
# with open(sf_all_pos) as fin_pos:
#     for line in fin_pos:
#         fields=line.split()
#         m_all_pos.append((fields[0], fields[1], fields[2]))
#
# m_support={}
# with open(sf_support_index) as fin_index:
#     for line in fin_index:
#         if line[0]==">":
#             fields=line.split("_")
#             m_support[int(fields[-1])]=1
#
#
# for i in range(len(m_all_pos)):
#     if m_support.has_key(i)==False:
#         print i, m_all_pos[i]

#Popen(cmd, shell = True, stdout = PIPE).communicate()

sf_1000g=sys.argv[1]
sf_repdenovo=sys.argv[2]

m_repdenovo={}
with open(sf_repdenovo,"r") as fin_repdenovo:
    for line in fin_repdenovo:
        if line[0]=="-":
            continue
        fields=line.split()
        pos=int(fields[0])
        m_repdenovo[pos]=1

cnt=0
all=0
with open(sf_1000g, "r") as fin_1000g:
    for line in fin_1000g:
        all=all+1
        fields=line.split()
        pos=int(fields[1])
        b_find=False
        for i in range(pos-20, pos+20):
            if m_repdenovo.has_key(i):
                cnt=cnt+1
                b_find=True
                break

        if b_find==True:
            print line.rstrip()

print cnt, all
