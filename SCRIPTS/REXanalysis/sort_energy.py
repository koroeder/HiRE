import numpy as np 
import os
import sys
import matplotlib.pyplot as plt

def main():
 
    T = np.loadtxt("temperatures.dat", dtype=float, usecols=0)
    T_list = list(np.round(T,4))
    nrep = len(T_list)
       
    fsize=np.loadtxt("md_energy.1.log")
    N = fsize.shape[0]
    U_kn = np.zeros([nrep, N])
    E_kn = np.zeros([nrep, N])
    T_kn = np.zeros([nrep, N])
    
    for i in range(1,nrep+1):
        U_kn[i-1]=np.loadtxt("md_energy."+str(i)+".log", dtype=float, usecols=2)
        T_kn[i-1]=np.loadtxt("md_temp."+str(i)+".log", dtype=float, usecols=1)
    
    for i in range(1,nrep+1):
        for k in range(N):
            posinT = T_list.index(T_kn[i-1][k])
            val = U_kn[posinT][k]
            # print(posinT, T_kn[i-1][k], val)
            E_kn[i-1][k] = val

      
    for i in range(1,nrep+1):
        np.savetxt("energy_sorted_temp."+str(i)+".log", E_kn[i-1])
#        np.savetxt("energy_sorted_rep."+str(i)+".log", U_kn[i-1])
        
        y, binEdges = np.histogram(E_kn[i-1][int(N/2):N], bins=50)
        bincenters = 0.5 * (binEdges[1:] + binEdges[:-1])
        plt.plot(bincenters, y, '-')
        plt.show()
    

if __name__ == "__main__":
  if len(sys.argv) > 1:
    main(int(sys.argv[1]))
  else:
    main()
