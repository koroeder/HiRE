import numpy as np
import sys

inp_x = sys.argv[1]
inp_pdb = sys.argv[2]
out_pdb = sys.argv[3]

coords = np.genfromtxt(inp_x, dtype = float)

outf = open(out_pdb,"w")

with open(inp_pdb,"r") as f:
    lines = f.readlines()
    for idx,line in enumerate(lines):
        newline = line[:30]
        newline += '{:8.3f}'.format(coords[idx,0])
        newline += '{:8.3f}'.format(coords[idx,1]) 
        newline += '{:8.3f}'.format(coords[idx,2])         
        newline += line[54:]
        outf.write(newline)
outf.close()

