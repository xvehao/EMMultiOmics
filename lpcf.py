#!/usr/bin/env python3.8
from __future__ import print_function
import sys
from pandas import read_csv
from mpmath import *
import os
#os.chdir("C:/Users/Yuanyuan Bian/Documents/2014 Spring/9530/final project")
mypath = sys.argv[1]
os.chdir(mypath)
a0 = sys.argv[2]
gstr = sys.argv[3]
D=read_csv("".join(["zv",str(a0),"_",gstr,".txt"]))
i=int(D.columns.tolist()[0])

result=[]
with open("".join(["zv",str(a0),"_",gstr,".txt"]),"r")as f:
    for line in f:
        result.append(list(map(float,line.split(","))))
        
v=result[0][1:(i+1)]
z=result[0][(i+1):(2*i+1)]

A=i*[0]

for j in range(i):
    A[j]=float(log(pcfd(v[j],z[j])))

out1=open("".join(["out1",str(a0),"_",gstr,".txt"]),"w")
print(" ".join( repr(e) for e in A ), file = out1)
out1.close()
