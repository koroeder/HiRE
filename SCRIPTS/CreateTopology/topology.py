# Data dictionary needed to create the topology
def init_topology():
    top_dict = dict()
    top_dict["nmol"] = 0        #number of molecules
    top_dict["nparticles"] = 0  #number of particles
    top_dict["ntypes"] = 0      #number of res types
    top_dict["npartypes"] = 0   #number of particle types
    top_dict["nbonds"] = 0      #number of bonds
    top_dict["nangles"] = 0     #number of angles
    top_dict["ndihs"] = 0       #number of dihedrals
    top_dict["nres"] = 0        #number of residues
    top_dict["nbondtypes"] = 0  #number of bond types
    top_dict["nangtypes"] = 0   #number of angular types
    top_dict["ndihstypes"] = 0  #number of dihedral types
    #now we have the different sections
    top_dict["particle_names"] = list()     #particle names
    top_dict["particle_mass"]  = list()     #particle mass
    top_dict["particle_type"] = list()      #particle type  
    top_dict["res_names"]  = list()         #residue names
    top_dict["res_start"] = list()          #residue first atom
    top_dict["res_final"] = list()          #residue final atom
    top_dict["rk"] = list()                 #bond spring constant
    top_dict["req"] = list()                #equilibrium bond length
    top_dict["tk"] = list()                 #angular spring constant
    top_dict["teq"]  = list()               #equilibrium angle
    top_dict["pk"] = list()                 #pk
    top_dict["pn"] = list()                 #pn
    top_dict["phi"]  = list()               #torsional phase
    top_dict["bonds"] = list()              #bond information
    top_dict["angs"]  = list()              #angle information
    top_dict["dihs"] = list()               #torsional information
    top_dict["fl_mol"] = list()             #chain pointers
    top_dict["charges"] = list()            #charges 
    return top_dict

# We use fix formatting to make parsing by fortran easier
# This function creates the fixed format entries 
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

#Function to create topology file from dictionary
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
    return