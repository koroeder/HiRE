## @file  ChemData.py
#
# @brief Contains functions to create chemical descriptors for topology
#
# Functions to obtain bond, angle and dihedral information
# The way these fucntions work is identical for all three in the algorithmic approach
# The functions for the bonding are annotated in detail

import NucleicAcidData

# @brief Function to obtain type from atom names
#
# Types are 1 indexed in the stored list of bonds in NucleicAcidData.
# If the type is used for the lookup later, -1 needs to be applied to each type index to get the correct zero-index.
def get_bondtype(at1,at2):
    for bond in NucleicAcidData.NA_bonds:
        if bond.iscorrectbond(at1, at2):
            return bond.get_typeid()
    return 0

# @brief Function to obtain angle type from atom names
def get_angletype(mol, at1, at2, at3):
    for angle in NucleicAcidData.NA_angles:
        if mol==angle.name[0:3]:
            if angle.iscorrectangle(at1, at2, at3):
                return angle.get_typeid()
    return 0

# @brief Function to obtain angle type from atom names
def get_dihtype(mol, at1, at2, at3, at4):
    types = list()
    for dih in NucleicAcidData.NA_dihs:
        if mol==dih.name[0:3]:
            if dih.iscorrectdih(at1, at2, at3, at4):
                types.append(dih.get_typeid())
                print(at1, at2, at3, at4)
    return types

# @brief Obtain bond information for topology
#
# Function obtains list of bonds, assuming it is a single molecule that is considered.\n
# We iterate over all particle names
# The bonds are fixed in our model, and hence we can find them by the first
# atom in the bond. It allows us to find all bonds and correctly identify 
# the parameters needed.\n
# bond_ids contains tuples of the bonded atom indices.\n
# bond_type contains the entry for the bond type found in RNA_bond from RNA_params.\n
# For each bond an entry is added to both lists.
#
# @param[in] CG_partnames - List of particle names for molecule

def get_bond_list(CG_partnames):
    bond_ids = list()
    bond_type = list()

    for idx,name in enumerate(CG_partnames):
        if name in ["P", "O", "C", "R4", "R1", "S1", "A1", "G1"]: #P-O5 bond, O5-C5 bond, C5-C4 bond, R4-R1 bond and R4-S1 bond, R1/S1-A/C/GU/T1 bond, A1-A2 and G1-G2 bond
            type = get_bondtype(name, CG_partnames[idx+1])
            if type > 0:
                bond_ids.append((idx,idx+1))
                bond_type.append(type - 1)
            else:
                raise ValueError("Error: Cannot find type for bond: ", name, CG_partnames[idx+1])
            if name == "R4":
                try:
                    if (CG_partnames[idx+3]) == "P": #R4-P bond
                        bond_ids.append((idx,idx+3))
                    elif (CG_partnames[idx+4]) == "P": #R4-P bond
                        bond_ids.append((idx,idx+4))
                    type = get_bondtype("R4", "P")
                    bond_type.append(type - 1)
                except IndexError:
                    continue
        else:
            continue
    return bond_ids, bond_type

## @brief get bonding information for all molecules
#
# To get the complete lists of bonds we go molecule by molecule.
# For each molecule, we call get_bond_list to get the bonds and type
# for the molecule. We record the offset for the first atom for each molecule,
# and then fix the atom ids if necessary.
# At the end, we have a full list of all bonds and types.
# Importantly, this information needs to be parsed into the correct
# format for the topology.
#
# @param[in] nmol - Number of molecules
# @param[in] termini - list of terminal atoms, used to derive the start and end indices for each molecule
# @param[in] CG_partnames - list of all CG particles names
#
# @see get_bond_list
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

## @brief Parsing bond information to different format
#
# Parsing the bond information produced so far into the correct format for 
# the topology file.\n
# The number of bond types is set (see RNA_params and DNA_params). In the topology, each bond type is refered to by an id,
# and the bond paramters are provided only once for each bond type.\n
# Here, the list of bonds and bondtypes is parsed into this format.
# 
# @param[in] bonds - list of bonded atoms as tuple
# @param[in] bondtype - list of bond types for each bond
# @param[in] moltype - list of type of molecule (currently RNA (0) or DNA (1)) to parse correct data
def get_bondinfo(bonds,bondtype):
    # Introduce bond types to track which ones we have already encountered
    bondtypes_used = [False for _ in range(len(NucleicAcidData.NA_bonds))] 
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
            #retrieve type of molecule
            bondtypes_used[this_type] = True
            nbondtypes += 1
            btypemap[this_type] = nbondtypes
            this_toptype = nbondtypes
            req.append(NucleicAcidData.NA_bonds[this_type].dist)
            rk.append(NucleicAcidData.NA_bonds[this_type].force)         
            bonds_top += [3*bond[0], 3*bond[1], this_toptype]           
    return nbondtypes,rk,req,nbonds,bonds_top


## @brief Obtain angle information for topology
#
# Function obtains list of bond angles, assuming it is a single molecule that is considered.\n
# We iterate over all particle names.
# The angles are fixed in our model, and hence we can find them by the first
# atom in the angle. It allows us to find all angles and correctly identify 
# the parameters needed.\n
# angle_ids contains tuples of the atom indices in the angle.\n
# angle_type contains the entry for the angle type found in RNA_angle from RNA_params.\n
# For each bond an entry is added to both lists.
#
# @param[in] CG_partnames - List of particle names for molecule
# @param[in] RNAorDNA - 0 - RNA, 1 - DNA
def get_angle_list(CG_partnames, RNAorDNA):
    angle_ids = list()
    angle_type = list()
    if RNAorDNA==0:
        moltype = "RNA"
    elif RNAorDNA==1:
        moltype = "DNA"
    for idx,name in enumerate(CG_partnames):
        if name in ["P", "O", "C", "R4"]: #P-O5-C5 angle, O5-C5-C4 angle, C5-C4-C1 angle, R4-R1-B1 angle 
            type = get_angletype(moltype, name, CG_partnames[idx+1], CG_partnames[idx+2])
            if type > 0:
                angle_ids.append((idx,idx+1,idx+2))
                angle_type.append(type - 1)             
            else:
                raise ValueError("Error: Cannot find type for bond: ", name, CG_partnames[idx+1])

            if name == "C": 
                try:
                    if (CG_partnames[idx+4]== "P"): #C5-C4-P angle
                        angle_ids.append((idx,idx+1,idx+4))
                        type = get_angletype(moltype,"C", "R4", "P")
                        angle_type.append(type - 1)
                    elif (CG_partnames[idx+5]== "P"):
                        angle_ids.append((idx,idx+1,idx+5))
                        type = get_angletype(moltype,"C", "R4", "P")
                        angle_type.append(type - 1)
                except IndexError:
                    continue

            elif name == "R4":
                try:
                    if CG_partnames[idx+3] == "P": #R4-P-O angle
                        angle_ids.append((idx,idx+3,idx+4))
                        type = get_angletype(moltype,"R4", "P", "O")
                        angle_type.append(type - 1)
                    elif CG_partnames[idx+4] == "P": #R4-P-O angle
                        angle_ids.append((idx,idx+4,idx+5))
                        type = get_angletype(moltype,"R4", "P", "O")
                        angle_type.append(type - 1)
                except IndexError:
                    continue
    
        elif name == "R1" or name == "S1":
            if CG_partnames[idx+1] == "A1": #R1-B1-B2 angle
                angle_ids.append((idx,idx+1,idx+2))
                type = get_angletype(moltype, name, CG_partnames[idx+1], CG_partnames[idx+2])
                angle_type.append(type - 1)
            elif CG_partnames[idx+1] == "G1":
                angle_ids.append((idx,idx+1,idx+2))
                type = get_angletype(moltype, name, CG_partnames[idx+1], CG_partnames[idx+2])
                angle_type.append(type - 1)
            try:
                if CG_partnames[idx+2] == "P": #R1-R4-P angle
                    angle_ids.append((idx,idx-1,idx+2))
                    type = get_angletype(moltype, "R1", "R4", "P")
                    angle_type.append(type - 1)
                elif CG_partnames[idx+3] == "P":
                    angle_ids.append((idx,idx-1,idx+3))
                    type = get_angletype(moltype, "R1", "R4", "P")
                    angle_type.append(type - 1)
            except IndexError:
                continue
        else:
            continue

    return angle_ids, angle_type

## @brief Parsing bond angle information to different format
#
# Parsing the bond angle information produced so far into the correct format for 
# the topology file.\n
# The number of bond angle types is set (see RNA_params and DNA_params). In the topology, each bond angle type is refered to by an id,
# and the bond angle paramters are provided only once for each bond angle type.\n
# Here, the list of angles and angle types is parsed into this format.
# 
# @param[in] angles - list of atoms in angle as tuple
# @param[in] angletype - list of bond angle types for each bond angle
# @param[in] moltype - list of type of molecule (currently RNA (0) or DNA (1)) to parse correct data
def get_angleinfo(angles,angletype):
    angletypes_used = [False for _ in range(len(NucleicAcidData.NA_angles))]
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
            #retrieve type of molecule
            angletypes_used[this_type] = True
            nangtypes += 1
            atypemap[this_type] = nangtypes
            this_toptype = nangtypes
            teq.append(NucleicAcidData.NA_angles[this_type].ang/180.0*3.1415926535897)
            tk.append(NucleicAcidData.NA_angles[this_type].force)         
            angles_top += [3*angle[0], 3*angle[1], 3*angle[2], this_toptype]           
    return nangtypes,tk,teq,nangles,angles_top



## @brief get bond angle information for all molecules
#
# To get the complete lists of bond angles we go molecule by molecule.
# For each molecule, we call get_angle_list to get the bond angles and type
# for the molecule. We record the offset for the first atom for each molecule,
# and then fix the atom ids if necessary.
# At the end, we have a full list of all bond angles and types.
# Importantly, this information needs to be parsed into the correct
# format for the topology.
#
# @param[in] nmol - Number of molecules
# @param[in] termini - list of terminal atoms, used to derive the start and end indices for each molecule
# @param[in] CG_partnames - list of all CG particles names
# @param[in] moltype - list of type of molecule (currently RNA (0) or DNA (1))
#
# @see get_angle_list
def get_angles(nmol, termini, CG_partnames, moltype):
    angles = list()
    angtype = list()
    for molid in range(nmol):
        offset = termini[2*molid] - 1
        start = termini[2*molid] - 1
        end = termini[2*molid+1]
        RNAorDNA = moltype[molid]
        angs_mol, angtype_mol = get_angle_list(CG_partnames[start:end],RNAorDNA)
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

## @brief Obtain dihedral information for topology
#
# Function obtains list of dihedrals, assuming it is a single molecule that is considered.\n
# We iterate over all particle names.
# The dihedrals are fixed in our model, and hence we can find them by the first
# atom. It allows us to find all dihedrals and correctly identify 
# the parameters needed.\n
# dih_ids contains tuples of the atom indices in the dihedral.\n
# dih_type contains the entry for the dihedral type found in RNA_dih from RNA_params.\n
# For each bond an entry is added to both lists.
#
# @param[in] CG_partnames - List of particle names for molecule
# @param[in] RNAorDNA - 0 - RNA, 1 - DNA
def get_dih_list(CG_partnames, RNAorDNA):
    dih_ids = list()
    dih_type = list()

    if RNAorDNA==0:
        moltype = "RNA"
    elif RNAorDNA==1:
        moltype = "DNA"
    for idx,name in enumerate(CG_partnames):
        if name in ["P", "O", "C"]: #P-O5-C5-C4 dih, O5-C5-C4-C1 dih, C5-C4-C1-B1 dih
            types = get_dihtype(moltype, name, CG_partnames[idx+1], CG_partnames[idx+2], CG_partnames[idx+3])
            for type in types:
                dih_ids.append((idx,idx+1,idx+2,idx+3))
                dih_type.append(type - 1)             

        if name == "P": 
            try:
                if ((CG_partnames[idx-2] == "A1") or (CG_partnames[idx-2] == "G1")): #P-R4-R1-B1 dih
                    types = get_dihtype(moltype, name, CG_partnames[idx-4], CG_partnames[idx-3], CG_partnames[idx-2])
                    for type in types:
                        dih_ids.append((idx,idx-4,idx-3,idx-2))
                        dih_type.append(type - 1)
                if ((CG_partnames[idx-1] == "C1") or (CG_partnames[idx-1] == "U1")): #P-R4-R1-B1 dih
                    types = get_dihtype(moltype, name, CG_partnames[idx-3], CG_partnames[idx-2], CG_partnames[idx-1])
                    for type in types:
                        dih_ids.append((idx,idx-3,idx-2,idx-1))
                        dih_type.append(type - 1)                    
            except IndexError:
                continue     

        elif name == "O": 
            try:
                if (CG_partnames[idx+5]== "P"): #O5-C5-C4-P dih
                    types = get_dihtype(moltype, name, CG_partnames[idx+1], CG_partnames[idx+2], CG_partnames[idx+5])
                    for type in types:
                        dih_ids.append((idx,idx+1,idx+2,idx+5))
                        dih_type.append(type - 1)
                elif (CG_partnames[idx+6]== "P"):
                    types = get_dihtype(moltype, name, CG_partnames[idx+1], CG_partnames[idx+2], CG_partnames[idx+6])
                    for type in types:
                        dih_ids.append((idx,idx+1,idx+2,idx+6))
                        dih_type.append(type - 1)                    
            except IndexError:
                continue

        elif name == "C":                        
            try:
                if (CG_partnames[idx+4]== "P"): #C5-C4-P-O dih
                    types = get_dihtype(moltype, name, CG_partnames[idx+1], CG_partnames[idx+4], CG_partnames[idx+5])
                    for type in types:
                        dih_ids.append((idx,idx+1,idx+4,idx+5))
                        dih_type.append(type - 1)                    
                elif (CG_partnames[idx+5]== "P"):
                    types = get_dihtype(moltype, name, CG_partnames[idx+1], CG_partnames[idx+5], CG_partnames[idx+6])
                    for type in types:
                        dih_ids.append((idx,idx+1,idx+5,idx+6))
                        dih_type.append(type - 1)                     
            except IndexError:
                continue

        elif name == "R4":          
            if ((CG_partnames[idx+2] == "A1") or (CG_partnames[idx+2] == "G1")): #R4-R1-B1-B2 dihs 
                types = get_dihtype(moltype, name, CG_partnames[idx+1], CG_partnames[idx+2], CG_partnames[idx+3])
                for type in types:
                    dih_ids.append((idx,idx+1,idx+2,idx+3))
                    dih_type.append(type - 1)
                types = get_dihtype(moltype, name, CG_partnames[idx+2], CG_partnames[idx+3], CG_partnames[idx+1])
                for type in types:
                    dih_ids.append((idx,idx+2,idx+3,idx+1))
                    dih_type.append(type - 1)      

            try:
                if CG_partnames[idx+3] == "P": #R4-P-O dih
                    types = get_dihtype(moltype, name, CG_partnames[idx+3], CG_partnames[idx+4], CG_partnames[idx+5])
                    for type in types:
                        dih_ids.append((idx,idx+3,idx+4,idx+5))
                        dih_type.append(type - 1)                    
                elif CG_partnames[idx+4] == "P": #R4-P-O dih
                    types = get_dihtype(moltype, name, CG_partnames[idx+4], CG_partnames[idx+5], CG_partnames[idx+6])
                    for type in types:
                        dih_ids.append((idx,idx+4,idx+5,idx+6))
                        dih_type.append(type - 1)                      
            except IndexError:
                continue

        elif name == "R1" or name == "S1":
            try:
                if CG_partnames[idx+2] == "P": #R1-R4-P-O dih
                    types = get_dihtype(moltype, name, CG_partnames[idx-1], CG_partnames[idx+2], CG_partnames[idx+3])
                    for type in types:
                        dih_ids.append((idx,idx-1,idx+2,idx+3))
                        dih_type.append(type - 1)
                elif CG_partnames[idx+3] == "P":
                    types = get_dihtype(moltype, name, CG_partnames[idx-1], CG_partnames[idx+3], CG_partnames[idx+4])
                    for type in types:
                        dih_ids.append((idx,idx-1,idx+3,idx+4))
                        dih_type.append(type - 1)                    
            except IndexError:
                continue           
        else:
            continue
     
    return dih_ids, dih_type



## @brief Parsing dihedral information to different format
#
# Parsing the dihedral information produced so far into the correct format for 
# the topology file.\n
# The number of dihedral types is set (see RNA_params and DNA_params). In the topology, each dihedral type is refered to by an id,
# and the dihedral paramters are provided only once for each dihedral type.\n
# Here, the list of dihedrals and dihedral types is parsed into this format.
# 
# @param[in] dihs - list of atoms in dihderals as tuples
# @param[in] dihtype - list of dihedral types for each dihedral
# @param[in] moltype - list of type of molecule (currently RNA (0) or DNA (1)) to parse correct data

def get_dihinfo(dihs,dihtype):
    dihtypes_used = [False for _ in range(len(NucleicAcidData.NA_dihs))]
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
            phi.append(NucleicAcidData.NA_dihs[this_type].dih/180.0*3.1415926535897)
            pk.append(NucleicAcidData.NA_dihs[this_type].force)
            pn.append(NucleicAcidData.NA_dihs[this_type].nterm)
            dihs_top += [3*dih[0], 3*dih[1], 3*dih[2], 3*dih[3], this_toptype]           
    return ntorstypes,pk,phi,ndihs,dihs_top,pn


## @brief get dihedral information for all molecules
#
# To get the complete lists of dihedrals we go molecule by molecule.
# For each molecule, we call get_dih_list to get the dihedrals and type
# for the molecule. We record the offset for the first atom for each molecule,
# and then fix the atom ids if necessary.
# At the end, we have a full list of all dihedrals and types.
# Importantly, this information needs to be parsed into the correct
# format for the topology.
#
# @param[in] nmol - Number of molecules
# @param[in] termini - list of terminal atoms, used to derive the start and end indices for each molecule
# @param[in] CG_partnames - list of all CG particles names
#
# @see get_dih_list
def get_dihs(nmol, termini, CG_partnames, moltype):
    dihs = list()
    torstype = list()
    for molid in range(nmol):
        offset = termini[2*molid] - 1
        start = termini[2*molid] - 1
        end = termini[2*molid+1]
        RNAorDNA = moltype[molid]
        angs_mol, torstype_mol = get_dih_list(CG_partnames[start:end],RNAorDNA)
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


def RNA_or_DNA(CG_resnames):
# check whether each molecule is RNA or DNA (can only be one!)
    nmol = 0
    termini_res = list()
    for idx,res in enumerate(CG_resnames):
        if res[-1] == "5":
            nmol += 1
            termini_res.append(idx)
        elif res[-1] == "3":
            termini_res.append(idx)           
    moltype = list()
    for mol in range(nmol):
        start = termini_res[2*mol]
        end = termini_res[2*mol+1]
        DNAT = False
        RNAT = False
       
        for res in CG_resnames[start:end+1]:
            if res in ["DA", "DA3", "DA5", "DC", "DC3", "DC5",
                       "DG", "DG3", "DG5", "DT", "DT3", "DT5"]:
                DNAT = True
            else:
                RNAT = True
        if (DNAT and RNAT):
            raise ValueError("This molecule has RNA and DNA nucleotides")
        elif (DNAT==0) and (RNAT==0):
            raise ValueError("This molecule has no RNA and DNA nucleotides")
        else:
            if RNAT:
                moltype.append(0)
            elif DNAT:
                moltype.append(1)
    return moltype