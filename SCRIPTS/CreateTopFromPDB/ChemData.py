# Functions to obtain bond, angle and dihedral information
# The way these fucntions work is identical for all three in the algorithmic approach
# The functions for the bonding are annotated in detail

import RNA_params
import DNA_params

# getting list of bonds, assuming it is a single molecule
def get_bond_list(CG_partnames):
    bond_ids = list()
    bond_type = list()
    # We iterate over all particle names
    # The bonds are fixed in our model, and hence we can find them by the first
    # atom in the bond. It allows us to find all bonds and correctly identify 
    # the parameters needed
    # bond_ids contains tuples of the bonded atom indices
    # bond_type contains the entry for the bond type found in RNA_bond from RNA_params
    # For each bond an entry is added to both lists
    for idx,name in enumerate(CG_partnames):
        if name == "P": #P-O5 bond
            bond_ids.append((idx,idx+1))
            bond_type.append(9)
        elif name == "O": #O5-C5 bond
            bond_ids.append((idx,idx+1))
            bond_type.append(10)
        elif name == "C": #C5-C4 bond
            bond_ids.append((idx,idx+1))
            bond_type.append(8)
        elif name == "R4":
            bond_ids.append((idx,idx+1)) #R4-R1 bond
            bond_type.append(1)
            try:
                if (CG_partnames[idx+3]) == "P": #R4-P bond
                    bond_ids.append((idx,idx+3))
                    bond_type.append(0)
                elif (CG_partnames[idx+4]) == "P": #R4-P bond
                    bond_ids.append((idx,idx+4))
                    bond_type.append(0)
            except IndexError:
                continue
        elif name == "R1":
            bond_ids.append((idx,idx+1)) #R1-B1 bond
            if CG_partnames[idx+1] == "A1":
                bond_type.append(3)
            elif CG_partnames[idx+1] == "C1":
                bond_type.append(5)
            elif CG_partnames[idx+1] == "G1":
                bond_type.append(2)
            elif CG_partnames[idx+1] == "U1":
                bond_type.append(4)       
        elif (name == "A1" or name == "G1"):
            bond_ids.append((idx,idx+1)) #B1-B2 bond
            if CG_partnames[idx+1] == "A2":
                bond_type.append(7)
            else:
                bond_type.append(6)
        else:
            continue
    return bond_ids, bond_type

# To get the complete lists of bonds we go molecule by molecule
# For each molecule, we call the above function to get the bonds and type
# for the molecule. We record the offset for the first atom for each molecule,
# and then fix the atom ids if necessary.
# At the end, we have a full list of all bonds and types.
# Importantly, this information needs to be parsed into the correct
# format for the topology.
def get_bonds(nmol, termini, CG_partnames):
    bonds = list()
    bondtype = list()
    for molid in range(nmol):
        offset = termini[2*molid] - 1
        start = termini[2*molid]
        end = termini[2*molid+1]
        bonds_mol, bondtype_mol = get_bond_list(CG_partnames[start-1:end])
        bondtype += bondtype_mol
        if offset == 0:
            bonds += bonds_mol           
        else:
            for bond in bonds_mol:
                at1 = bond[0]
                at2 = bond[1]
                bonds.append((at1+offset, at2+offset))
    return bonds,bondtype

# Parsing the bond information produced so far into the correct format for 
# the topology file. 
def get_bondinfo(bonds,bondtype):
    # Introduce bond types to track which ones we have already encountered
    bondtypes_used = [False, False, False, False, False, False, 
                      False, False, False, False, False]
    nbonds = len(bondtype)
    btypemap = dict()
    req = list()
    rk = list()
    bonds_top = list()
    nbondtypes = 0
    # Iterate over all the bonds found earlier
    for idx,bond in enumerate(bonds):
        # Retrieve the bond type
        this_type = bondtype[idx]
        # If this type has been found before, we simply add
        # the bond information (atom ids + type for topology)
        if bondtypes_used[this_type]:
            this_toptype = btypemap[this_type]
            bonds_top += [3*bond[0], 3*bond[1], this_toptype]
        # Otherwise, we add entries to the bond type maps,
        # entries to the spring constants and eq distance and 
        # change the number of bond types
        else:
            bondtypes_used[this_type] = True
            nbondtypes += 1
            btypemap[this_type] = nbondtypes
            this_toptype = nbondtypes
            req.append(RNA_params.RNA_bonds[this_type].dist)
            rk.append(RNA_params.RNA_bonds[this_type].force)
            bonds_top += [3*bond[0], 3*bond[1], this_toptype]           
    return nbondtypes,rk,req,nbonds,bonds_top


#getting list of angles, assuming it is a single molecule
def get_angle_list(CG_partnames):
    angle_ids = list()
    angle_type = list()
    for idx,name in enumerate(CG_partnames):
        if name == "P": #P-O5-C5 angle
            angle_ids.append((idx,idx+1,idx+2))
            angle_type.append(6)
        elif name == "O": #O5-C5-C4 angle
            angle_ids.append((idx,idx+1,idx+2))
            angle_type.append(7)
        elif name == "C": 
            angle_ids.append((idx,idx+1,idx+2)) #C5-C4-C1 angle
            angle_type.append(10)
            try:
                if (CG_partnames[idx+4]== "P"): #C5-C4-P angle
                    angle_ids.append((idx,idx+1,idx+4))
                    angle_type.append(8)
                elif (CG_partnames[idx+5]== "P"):
                    angle_ids.append((idx,idx+1,idx+5))
                    angle_type.append(8)
            except IndexError:
                continue
        elif name == "R4":
            angle_ids.append((idx,idx+1,idx+2)) #R4-R1-B1 bond           
            if CG_partnames[idx+2] == "A1":
                angle_type.append(0)
            elif CG_partnames[idx+2] == "C1":
                angle_type.append(3)
            elif CG_partnames[idx+2] == "G1":
                angle_type.append(2)
            elif CG_partnames[idx+2] == "U1":
                angle_type.append(1) 
            try:
                if CG_partnames[idx+3] == "P": #R4-P-O angle
                    angle_ids.append((idx,idx+3,idx+4))
                    angle_type.append(9)
                elif CG_partnames[idx+4] == "P": #R4-P-O angle
                    angle_ids.append((idx,idx+4,idx+5))
                    angle_type.append(9)
            except IndexError:
                continue
        elif name == "R1":
            if CG_partnames[idx+1] == "A1": #R1-B1-B2 angle
                angle_ids.append((idx,idx+1,idx+2))
                angle_type.append(4)
            elif CG_partnames[idx+1] == "G1":
                angle_ids.append((idx,idx+1,idx+2))
                angle_type.append(5)
            try:
                if CG_partnames[idx+2] == "P": #R1-R4-P angle
                    angle_ids.append((idx,idx-1,idx+2))
                    angle_type.append(11)
                elif CG_partnames[idx+3] == "P":
                    angle_ids.append((idx,idx-1,idx+3))
                    angle_type.append(11)
            except IndexError:
                continue
            
        else:
            continue
    return angle_ids, angle_type


def get_angleinfo(angles,angletype):
    angletypes_used = [False, False, False, False, False, False, 
                       False, False, False, False, False, False]
    nangles = len(angletype)
    atypemap = dict()
    teq = list()
    tk = list()
    angles_top = list()
    nangtypes = 0
    for idx,angle in enumerate(angles):
        this_type = angletype[idx]
        if angletypes_used[this_type]:
            this_toptype = atypemap[this_type]
            angles_top += [3*angle[0], 3*angle[1], 3*angle[2], this_toptype]
        else:
            angletypes_used[this_type] = True
            nangtypes += 1
            atypemap[this_type] = nangtypes
            this_toptype = nangtypes
            teq.append(RNA_params.RNA_angles[this_type].ang/180.0*3.1415926535897)
            tk.append(RNA_params.RNA_angles[this_type].force)
            angles_top += [3*angle[0], 3*angle[1], 3*angle[2], this_toptype]           
    return nangtypes,tk,teq,nangles,angles_top

def get_angles(nmol, termini, CG_partnames):
    angles = list()
    angtype = list()
    for molid in range(nmol):
        offset = termini[2*molid] - 1
        start = termini[2*molid] - 1
        end = termini[2*molid+1]
        angs_mol, angtype_mol = get_angle_list(CG_partnames[start:end])
        angtype += angtype_mol
        if offset == 0:
            angles += angs_mol           
        else:
            for ang in angs_mol:
                at1 = ang[0]
                at2 = ang[1]
                at3 = ang[2]
                angles.append((at1+offset, at2+offset, at3+offset))
    return angles,angtype

#getting list of dihs, assuming it is a single molecule
def get_dih_list(CG_partnames):
    dih_ids = list()
    dih_type = list()
    for idx,name in enumerate(CG_partnames):
        if name == "P": #P-O5-C5-C4 dih
            dih_ids.append((idx,idx+1,idx+2,idx+3))
            dih_type.append(21)
            dih_ids.append((idx,idx+1,idx+2,idx+3))
            dih_type.append(22) 
            dih_ids.append((idx,idx+1,idx+2,idx+3))
            dih_type.append(23)  
            if CG_partnames[idx-2] == "A1": #P-R4-R1-B1 dih
                dih_ids.append((idx,idx-4,idx-3,idx-2))
                dih_type.append(13)
            if CG_partnames[idx-1] == "C1":
                dih_ids.append((idx,idx-3,idx-2,idx-1))
                dih_type.append(15)  
            if CG_partnames[idx-2] == "G1":
                dih_ids.append((idx,idx-4,idx-3,idx-2))
                dih_type.append(14)
            if CG_partnames[idx-1] == "U1":
                dih_ids.append((idx,idx-3,idx-2,idx-1))
                dih_type.append(16)                 
        elif name == "O": #O5-C5-C4-C1 dih
            dih_ids.append((idx,idx+1,idx+2,idx+3))
            dih_type.append(20)
            try:
                if (CG_partnames[idx+5]== "P"): #O5-C5-C4-P dih
                    dih_ids.append((idx,idx+1,idx+2,idx+5))
                    dih_type.append(19)
                elif (CG_partnames[idx+6]== "P"):
                    dih_ids.append((idx,idx+1,idx+2,idx+6))
                    dih_type.append(19)
            except IndexError:
                continue
        elif name == "C":
            #C5-C4-C1-B1 dih
            if (CG_partnames[idx+3] == "A1"):
                dih_ids.append((idx,idx+1,idx+2,idx+3)) 
                dih_ids.append((idx,idx+1,idx+2,idx+3))
                dih_type.append(6)
                dih_type.append(7)
            elif (CG_partnames[idx+3] == "C1"):
                dih_ids.append((idx,idx+1,idx+2,idx+3)) 
                dih_ids.append((idx,idx+1,idx+2,idx+3))
                dih_type.append(9)
                dih_type.append(10)   
            elif (CG_partnames[idx+3] == "G1"):
                dih_ids.append((idx,idx+1,idx+2,idx+3)) 
                dih_type.append(8)
            elif (CG_partnames[idx+3] == "U1"):
                dih_ids.append((idx,idx+1,idx+2,idx+3)) 
                dih_ids.append((idx,idx+1,idx+2,idx+3))
                dih_type.append(11)
                dih_type.append(12)                               
            try:
                if (CG_partnames[idx+4]== "P"): #C5-C4-P-O dih
                    dih_ids.append((idx,idx+1,idx+4,idx+5))
                    dih_type.append(17)
                elif (CG_partnames[idx+5]== "P"):
                    dih_ids.append((idx,idx+1,idx+5,idx+6))
                    dih_type.append(17)
            except IndexError:
                continue
        elif name == "R4":          
            if CG_partnames[idx+2] == "A1": #R4-R1-B1-B2 dihs 
                dih_ids.append((idx,idx+1,idx+2,idx+3)) 
                dih_type.append(2)
                dih_ids.append((idx,idx+1,idx+2,idx+3)) 
                dih_type.append(3)
                dih_ids.append((idx,idx+2,idx+3,idx+1))
                dih_type.append(4)
            elif CG_partnames[idx+2] == "G1": #R4-R1-B1-B2 dihs 
                dih_ids.append((idx,idx+1,idx+2,idx+3)) 
                dih_type.append(0)
                dih_ids.append((idx,idx+1,idx+2,idx+3)) 
                dih_type.append(1)
                dih_ids.append((idx,idx+2,idx+3,idx+1))
                dih_type.append(5)
            try:
                if CG_partnames[idx+3] == "P": #R4-P-O dih
                    dih_ids.append((idx,idx+3,idx+4,idx+5))
                    dih_type.append(24)
                elif CG_partnames[idx+4] == "P": #R4-P-O dih
                    dih_ids.append((idx,idx+4,idx+5,idx+6))
                    dih_type.append(24)
            except IndexError:
                continue
        elif name == "R1":
            try:
                if CG_partnames[idx+2] == "P": #R1-R4-P-O dih
                    dih_ids.append((idx,idx-1,idx+2,idx+3))
                    dih_type.append(18)
                elif CG_partnames[idx+3] == "P":
                    dih_ids.append((idx,idx-1,idx+3,idx+4))
                    dih_type.append(18)
            except IndexError:
                continue           
        else:
            continue
    return dih_ids, dih_type



def get_dihinfo(dihs,dihtype):
    dihtypes_used = [False, False, False, False, False,
                     False, False, False, False, False, 
                     False, False, False, False, False,
                     False, False, False, False, False,
                     False, False, False, False, False]
    ndihs = len(dihtype)
    atypemap = dict()
    phi = list()
    pk = list()
    dihs_top = list()
    pn = list()
    ntorstypes = 0
    for idx,dih in enumerate(dihs):
        this_type = dihtype[idx]
        if dihtypes_used[this_type]:
            this_toptype = atypemap[this_type]
            dihs_top += [3*dih[0], 3*dih[1], 3*dih[2], 3*dih[3], this_toptype]
        else:
            dihtypes_used[this_type] = True
            ntorstypes += 1
            atypemap[this_type] = ntorstypes
            this_toptype = ntorstypes
            phi.append(RNA_params.RNA_dihs[this_type].dih/180.0*3.1415926535897)
            pk.append(RNA_params.RNA_dihs[this_type].force)
            pn.append(RNA_params.RNA_dihs[this_type].nterm)
            dihs_top += [3*dih[0], 3*dih[1], 3*dih[2], 3*dih[3], this_toptype]           
    return ntorstypes,pk,phi,ndihs,dihs_top,pn

def get_dihs(nmol, termini, CG_partnames):
    dihs = list()
    torstype = list()
    for molid in range(nmol):
        offset = termini[2*molid] - 1
        start = termini[2*molid] - 1
        end = termini[2*molid+1]
        angs_mol, torstype_mol = get_dih_list(CG_partnames[start:end])
        torstype += torstype_mol
        if offset == 0:
            dihs += angs_mol           
        else:
            for ang in angs_mol:
                at1 = ang[0]
                at2 = ang[1]
                at3 = ang[2]
                at4 = ang[3]
                dihs.append((at1+offset, at2+offset, at3+offset, at4+offset))
    
    return dihs,torstype
