class particle:
    def __init__(self, name, el, idx, charge, mass):
        self.name = name
        self.el = el
        self.idx = idx
        self.charge = charge
        self.mass = mass
    
    def set_el(self,newel):
        self.el = newel
    
    def set_charge(self,newq):
        self.charge = newq
    
    def set_name(self,newname):
        self.name = newname
    
    def set_mass(self,newmass):
        self.mass = newmass
    
class residue:
    def __init__(self, rtype, idx, start, end):
        self.rtype = rtype
        self.idx = idx
        self.start = start
        self.end = end
        self.npar = end - start + 1

class molecule(particle,residue):
    def __init__(self, nres, npars, pfirst, plast, rfirst, rlast, idx):
        self.nres = nres       #number of residues
        self.npars = npars     #number of particles
        self.pfirst = pfirst   #first particle
        self.plast = plast     #last particle
        self.rfirst = rfirst   #first residue
        self.rlast = rlast     #last residue
        self.idx = idx         #index
        self.reslist = list()
        self.atomlist = list()
    
    def set_nres(self,newnres):
        self.nres = newnres
    

class system(molecule,particle,residue):
    def __init__(self,nmol,nres):
        self.nmol = nmol
        self.nres = nres
        self.mollist = list()
        self.reslist = list()
        
    def create_allmols(self, first, last, rfirst, rlast):
        for idx in range(self.nmol):
            nparticles = last[idx] - first[idx] + 1
            nres = rlast[idx] - rfirst[idx] + 1
            self.mollist.append(molecule(nres, nparticles, first[idx], 
                                         last[idx], rfirst[idx], rlast[idx], 
                                         idx+1))
    
    def add_reslist(self,restype,first,last):
        for idx,(r,f,l) in enumerate(zip(restype,first,last)):
            self.reslist.append(residue(r,idx,f,l))

    
    
    