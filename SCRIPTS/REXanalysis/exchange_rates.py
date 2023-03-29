import numpy as np 
import os
import sys
import matplotlib.pyplot as plt

filename = sys.argv[1]
nrep = int(sys.argv[2])

rates = np.zeros(nrep)
excount = 0

outfile = open(filename, "r")
repline = outfile.readlines()

for i in range(len(repline)):
    if "sel_exchanges" in repline[i] :
        excount += 1
    
    if "T_REX> exchanged replicas" in repline[i] : 
        ex1 = int(repline[i].split()[3])
        ex2 = int(repline[i].split()[5])
        #print(ex1, ex2)
        rates[ex1-1] +=1
        rates[ex2-1] +=1

print("echagnges", excount)
for i in rates:
    print('{:.3f}'.format(i/excount))
                  
        
