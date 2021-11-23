import topology_classes as top
import math

#input is filename for ichain, 
#output is number of molecules and their first and last particle as a list
def parse_ichain(fname):
    #open ichain_RNA.dat file
    inunit = open(fname, "r")
    #first line is the number of molecules
    nmols = int(inunit.readline().split()[0])
    #read information about each molecule
    m_pfirst = list()
    m_pfinal = list()
    for l_id in range(nmols):
        l = inunit.readline().split()
        m_pfirst.append(int(l[0]))
        m_pfinal.append(int(l[1]))
    return nmols,m_pfirst,m_pfinal


#parse baselist.dat
def parse_baselist(fname):
    r_pfinal = list()
    restype = list()
    nres = 0
    #open baselist.dat
    with open(fname, "r") as f:
        for line in f:
            nres += 1
            l = line.split()
            r_pfinal.append(int(l[0]))
            restype.append(int(l[1]))
    return nres,r_pfinal,restype

#get base information (first and last particle and type)
def get_baseinfo(fname):
    nres,r_pfinal,restype = parse_baselist(fname)
    r_pfirst = list()
    for rid in range(nres):
        if rid == 0:
            r_pfirst.append(1)
        else:
            r_pfirst.append(r_pfinal[rid-1]+1)
    return nres,r_pfirst,r_pfinal,restype

#iterate over the list of residues and first and last particles in the 
#molecules to get the first and last residue in each molecule
def create_respointers_mol(nmols,m_pfirst,m_pfinal,r_pfirst,r_pfinal):
    m_rfirst = list()
    m_rfinal = list()
    for m_id in range(nmols):
        start = m_pfirst[m_id]
        end = m_pfinal[m_id]
        for idx,(s,e) in enumerate(zip(r_pfirst,r_pfinal)):
            if start==s:
                m_rfirst.append(idx+1)
            if end==e:
                m_rfinal.append(idx+1)
    return m_rfirst,m_rfinal


def get_mol_base_info(fchain="ichain_RNA.dat",fbases="baselist.dat"):
    nmols,m_pfirst,m_pfinal = parse_ichain(fchain)
    nres,r_pfirst,r_pfinal,restype = get_baseinfo(fbases)
    m_rfirst,m_rfinal = create_respointers_mol(nmols, m_pfirst, m_pfinal, 
                                               r_pfirst, r_pfinal)
    systemRNA = top.system(nmols,nres)
    systemRNA.create_allmols(m_pfirst, m_pfinal, m_rfirst, m_rfinal)
    systemRNA.add_reslist(restype, r_pfirst, r_pfinal)
    return systemRNA

def parse_chargedat(fname):
    pcharge = list()
    with open(fname, "r") as f:
        for line in f:
            l = line.split()
            pcharge.append(int(l[1]))
    return pcharge

#read section assuming entries are space separated
def read_section(nentries,entriesperline,unit):
    nlines = math.ceil(float(nentries)/entriesperline)   
    out = list()
    for i in range(nlines):
        line = unit.readline().replace("\n","").replace("\r","")     
        out += line.split()
    return out

#read fixed width section,where entries may be fused
def read_section_fw(nentries,entriesperline,entrywidth,unit):
    nlines = math.ceil(float(nentries)/entriesperline)
    out = list()
    for i in range(nlines):
        line = unit.readline().replace("\n","").replace("\r","")
        list_of_strings = list()
        for i in range(0, len(line), entrywidth):
            list_of_strings.append(line[i:entrywidth+i].strip())
        out += list_of_strings    
    return out
 
def parse_oldtop(fname):
    top_data = dict()
    topunit = open(fname, "r")
    #first line is the name, we simply skip ahead
    topunit.readline()
    #the second line contains the sizes of arrays etc.
    line = topunit.readline().split()
    npar = line[0]              #number of particles
    ntypes = line[1]            #number of res types
    npartypes = line[2]         #number of particle types
    nbonds = line[3]            #number of bonds
    nangles = line[4]           #number of angles
    ndihs = line[5]             #number of dihedrals
    nres = line[6]              #number of residues
    nbondtypes = line[7]        #number of bond types
    nangtypes = line[8]         #number of angular types
    ndihstypes = line[9]        #number of dihedral types
    top_data["nparticles"] = npar
    top_data["ntypes"] = ntypes
    top_data["npartypes"] = npartypes
    top_data["nbonds"] = nbonds
    top_data["nangles"] = nangles
    top_data["ndihs"] = ndihs
    top_data["nres"] = nres
    top_data["nbondtypes"] = nbondtypes
    top_data["nangtypes"] =  nangtypes
    top_data["ndihstypes"] =  ndihstypes
    #now we have the different sections
    #particle names
    top_data["particle_names"] = read_section(int(npar),20,topunit)
    #particle mass
    top_data["particle_mass"]= [float(m) for m in read_section(int(npar),5,topunit)]  
    #particle type
    top_data["particle_type"]= [int(t) for t in read_section(int(npar),12,topunit)]  
    #residue names
    top_data["res_names"] = read_section_fw(int(nres),20,4,topunit)  
    #residue first atom
    top_data["res_start"] = [int(s) for s in read_section(int(nres),12,topunit)] 
    #bond spring constant
    top_data["rk"] = [float(rk) for rk in read_section(int(nbondtypes),5,topunit)]
    #equilibrium bond length
    top_data["req"] = [float(re) for re in read_section(int(nbondtypes),5,topunit)]
    #angular spring constant
    top_data["tk"] = [float(tk) for tk in read_section(int(nangtypes),5,topunit)]
    #equilibrium angle
    top_data["teq"] = [float(te) for te in read_section(int(nangtypes),5,topunit)]  
    #pk
    top_data["pk"] = [float(pk) for pk in read_section(int(ndihstypes),5,topunit)]
    #pn
    top_data["pn"] = [float(pn) for pn in read_section(int(ndihstypes),5,topunit)]
    #torsional phase
    top_data["phi"]= [float(phi) for phi in read_section(int(ndihstypes),5,topunit)]   
    #bond information
    top_data["bonds"] = [int(i) for i in read_section(3*int(nbonds),12,topunit)]
    #angle information
    top_data["angs"] = [int(i) for i in read_section(4*int(nangles),12,topunit)]  
    #torsional information
    top_data["dihs"] = [int(i) for i in read_section(5*int(ndihs),12,topunit)]
    topunit.close()
    return top_data


def get_top_dict(systemRNA,fname="parametres_RNA.top"):
    top_dict = parse_oldtop(fname)
    top_dict["nmol"] = systemRNA.nmol
    #now get first and last particles for each fragment
    pfl = list()
    for mol in systemRNA.mollist:
        pfl.append(mol.pfirst)
        pfl.append(mol.plast)
    top_dict["fl_mol"] = pfl
    #get res_final
    res_end = list()
    for res in systemRNA.reslist:
        res_end.append(res.end)
    top_dict["res_final"] = res_end
    top_dict["charges"] = parse_chargedat("chargeatm_RNA.dat")
    return top_dict

def create_fw_lines(data,nentries,formatstr):
    lines = list()
    line = ""
    for idx,d in enumerate(data):
        if ((idx+1)%nentries)!=0:
            line += formatstr.format(d)
        else:
            line += formatstr.format(d)
            line += "\n"
            lines.append(line)
            line = ""
    if line!="":
        line += "\n"
        lines.append(line)
    return lines

def create_newtopology(top_dict,fname="parameters.top"):
    outf = open(fname, "w")
    outf.write("Molecule RNA\n")
    #write definitions
    outf.write("SECTION DEFINITIONS\n")
    data = [top_dict["nparticles"], top_dict["npartypes"], top_dict["nbonds"],
            top_dict["nangles"], top_dict["ndihs"], top_dict["nres"], 
            top_dict["nbondtypes"], top_dict["nangtypes"], 
            top_dict["ndihstypes"], str(top_dict["nmol"])]
    lines = create_fw_lines(data, 12,'{:>6}')
    for line in lines:
        outf.write(line)
    #write particle names
    outf.write("SECTION PARTICLE_NAMES\n")
    data = top_dict["particle_names"]
    lines = create_fw_lines(data, 20,'{:4}')
    for line in lines:
        outf.write(line)    
    #write residue labels
    outf.write("SECTION RESIDUE_LABELS\n")
    data = top_dict["res_names"]
    lines = create_fw_lines(data, 20,'{:4}')
    for line in lines:
        outf.write(line)  
    #write residue start and finish ids
    outf.write("SECTION RESIDUE_POINTER\n")
    data = list()
    for s,e in zip(top_dict["res_start"],top_dict["res_final"]):
        data += [s,e]
    lines = create_fw_lines(data, 12,'{:6d}')
    for line in lines:
        outf.write(line)         
    #write chain pointer
    outf.write("SECTION CHAIN_POINTER\n")
    data = top_dict["fl_mol"]
    lines = create_fw_lines(data, 12,'{:6d}')
    for line in lines:
        outf.write(line)   
    #write particle mass
    outf.write("SECTION PARTICLE_MASSES\n")
    data = top_dict["particle_mass"]
    lines = create_fw_lines(data, 5,'{:16.8f}')
    for line in lines:
        outf.write(line)          
    #write 
    outf.write("SECTION PARTICLE_TYPE\n")
    data = top_dict["particle_type"]
    lines = create_fw_lines(data, 12,'{:6d}')
    for line in lines:
        outf.write(line) 
    #write
    outf.write("SECTION CHARGES\n")
    data = top_dict["charges"]
    lines = create_fw_lines(data, 5,'{:16.8f}')
    for line in lines:
        outf.write(line)   
    #write
    outf.write("SECTION BOND_FORCE_CONSTANT\n")
    data = top_dict["rk"]
    lines = create_fw_lines(data, 5,'{:16.8f}')
    for line in lines:
        outf.write(line)  
    #write
    outf.write("SECTION BOND_EQUIL_VALUE\n")
    data = top_dict["req"]
    lines = create_fw_lines(data, 5,'{:16.8f}')
    for line in lines:
        outf.write(line)  
    #write
    outf.write("SECTION ANGLE_FORCE_CONSTANT\n")
    data = top_dict["tk"]
    lines = create_fw_lines(data, 5,'{:16.8f}')
    for line in lines:
        outf.write(line)          
    #write
    outf.write("SECTION ANGLE_EQUIL_VALUE\n")
    data = top_dict["teq"]
    lines = create_fw_lines(data, 5,'{:16.8f}')
    for line in lines:
        outf.write(line)          
    #write
    outf.write("SECTION DIHEDRAL_FORCE_CONSTANT\n")
    data = top_dict["pk"]
    lines = create_fw_lines(data, 5,'{:16.8f}')
    for line in lines:
        outf.write(line)  
    #write
    outf.write("SECTION DIHEDRAL_PERIODICITY\n")
    data = top_dict["pn"]
    lines = create_fw_lines(data, 5,'{:16.8f}')
    for line in lines:
        outf.write(line)  
    #write
    outf.write("SECTION DIHEDRAL_PHASE\n")
    data = top_dict["phi"]
    lines = create_fw_lines(data, 5,'{:16.8f}')
    for line in lines:
        outf.write(line)  
    #write 
    outf.write("SECTION BONDS\n")
    data = top_dict["bonds"]
    lines = create_fw_lines(data, 12,'{:6d}')
    for line in lines:
        outf.write(line)         
    #write 
    outf.write("SECTION ANGLES\n")
    data = top_dict["angs"]
    lines = create_fw_lines(data, 12,'{:6d}')
    for line in lines:
        outf.write(line)         
    #write 
    outf.write("SECTION DIHEDRALS\n")
    data = top_dict["dihs"]
    lines = create_fw_lines(data, 12,'{:6d}')
    for line in lines:
        outf.write(line) 

    outf.close()


if __name__ == "__main__":
    systemRNA = get_mol_base_info()
    top_dict = get_top_dict(systemRNA)
    create_newtopology(top_dict)
    
    
