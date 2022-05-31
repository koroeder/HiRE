import RNA_params
import DNA_params

#getting list of bonds, assuming it is a single molecule
def get_bond_list(CG_partnames):
    bond_ids = list()
    bond_type = list()
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


def get_bondinfo(bonds,bondtype):
    bondtypes_used = [False, False, False, False, False, False, 
                      False, False, False, False, False]
    nbonds = len(bondtype)
    btypemap = dict()
    req = list()
    rk = list()
    bonds_top = list()
    nbondtypes = 0
    for idx,bond in enumerate(bonds):
        this_type = bondtype[idx]
        if bondtypes_used[this_type]:
            this_toptype = btypemap[this_type]
            bonds_top += [3*bond[0], 3*bond[1], this_toptype]
        else:
            bondtypes_used[this_type] = True
            nbondtypes += 1
            btypemap[this_type] = nbondtypes
            this_toptype = nbondtypes
            req.append(RNA_params.RNA_bonds[this_type].dist)
            rk.append(RNA_params.RNA_bonds[this_type].force)
            bonds_top += [3*bond[0], 3*bond[1], this_toptype]           
    return nbondtypes,rk,req,nbonds,bonds_top

def get_bonds(nmol, termini, CG_partnames):
    bonds = list()
    bondtype = list()
    for molid in range(nmol):
        offset = termini[2*molid] - 1
        start = termini[2*molid]
        end = termini[2*molid+1] + 1
        bonds_mol, bondtype_mol = get_bond_list(CG_partnames[start:end])
        bondtype += bondtype_mol
        if offset == 0:
            bonds += bonds_mol           
        else:
            for bond in bonds_mol:
                at1 = bond[0]
                at2 = bond[1]
                bonds.append((at1+offset, at2+offset))
    
    return bonds,bondtype


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
            angle_type.append(1)            
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
            teq.append(RNA_params.RNA_angles[this_type].ang)
            tk.append(RNA_params.RNA_angles[this_type].force)
            angles_top += [3*angle[0], 3*angle[1], 3*angle[2], this_toptype]           
    return nangtypes,tk,teq,nangles,angles_top

def get_angles(nmol, termini, CG_partnames):
    angles = list()
    angtype = list()
    for molid in range(nmol):
        offset = termini[2*molid] - 1
        start = termini[2*molid]
        end = termini[2*molid+1] + 1
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

    #continue here
    6	Dihedral_RNA( "C", "R4", "R1", "A1", 1.0,  -150.0, 1),
	7	Dihedral_RNA( "C", "R4", "R1", "A1", 0.2,   150.0, 1),
	8	Dihedral_RNA( "C", "R4", "R1", "G1", 1.0,  -150.0, 1),
	9	Dihedral_RNA( "C", "R4", "R1", "C1", 1.0,  -150.0, 1),
	10	Dihedral_RNA( "C", "R4", "R1", "C1", 0.2,   150.0, 1),
	11	Dihedral_RNA( "C", "R4", "R1", "U1", 1.0,  -150.0, 1),
	12	Dihedral_RNA( "C", "R4", "R1", "U1", 0.2,   150.0, 1),
	17	Dihedral_RNA( "C", "R4",  "P",  "O", 1.2,    20.0, 1),
        elif name == "C": 
            dih_ids.append((idx,idx+1,idx+2)) #C5-C4-C1 dih
            dih_type.append(10)
            try:
                if (CG_partnames[idx+4]== "P"): #C5-C4-P dih
                    dih_ids.append((idx,idx+1,idx+4))
                    dih_type.append(8)
                elif (CG_partnames[idx+5]== "P"):
                    dih_ids.append((idx,idx+1,idx+5))
                    dih_type.append(8)
            except IndexError:
                continue
    
    0   Dihedral_RNA("R4", "R1", "G1", "G2", 1.0,     -20, 1),
	1	Dihedral_RNA("R4", "R1", "G1", "G2", 0.2,    -160, 1),
	2	Dihedral_RNA("R4", "R1", "A1", "A2", 1.0,     -20, 1),
	3	Dihedral_RNA("R4", "R1", "A1", "A2", 0.2,  -160.0, 1),    
	4	Dihedral_RNA("R4", "A1", "A2", "R1", 1.0,  -165.0, 1),
	5	Dihedral_RNA("R4", "G1", "G2", "R1", 1.0,  -165.0, 1),
	24	Dihedral_RNA("R4",  "P",  "O",  "C", 1.00,    0.0, 2)] 
        elif name == "R4":
            dih_ids.append((idx,idx+1,idx+2)) #R4-R1-B1 bond
            dih_type.append(1)            
            if CG_partnames[idx+2] == "A1":
                dih_type.append(0)
            elif CG_partnames[idx+2] == "C1":
                dih_type.append(3)
            elif CG_partnames[idx+2] == "G1":
                dih_type.append(2)
            elif CG_partnames[idx+2] == "U1":
                dih_type.append(1) 
            try:
                if CG_partnames[idx+3] == "P": #R4-P-O dih
                    dih_ids.append((idx,idx+3,idx+4))
                    dih_type.append(9)
                elif CG_partnames[idx+4] == "P": #R4-P-O dih
                    dih_ids.append((idx,idx+4,idx+5))
                    dih_type.append(9)
            except IndexError:
                continue

    	18	Dihedral_RNA("R1", "R4",  "P",  "O", 1.2,   160.0, 1),
        elif name == "R1":
            if CG_partnames[idx+1] == "A1": #R1-B1-B2 dih
                dih_ids.append((idx,idx+1,idx+2))
                dih_type.append(4)
            elif CG_partnames[idx+1] == "G1":
                dih_ids.append((idx,idx+1,idx+2))
                dih_type.append(5)
            try:
                if CG_partnames[idx+2] == "P": #R1-R4-P dih
                    dih_ids.append((idx,idx-1,idx+2))
                    dih_type.append(11)
                elif CG_partnames[idx+3] == "P":
                    dih_ids.append((idx,idx-1,idx+3))
                    dih_type.append(11)
            except IndexError:
                continue
            
        else:
            continue
    return dih_ids, dih_type



def get_dihinfo(dihs,dihtype):
    dihtypes_used = [False, False, False, False, False, False, 
                       False, False, False, False, False, False]
    ndihs = len(dihtype)
    atypemap = dict()
    teq = list()
    tk = list()
    dihs_top = list()
    ntorstypes = 0
    for idx,dih in enumerate(dihs):
        this_type = dihtype[idx]
        if dihtypes_used[this_type]:
            this_toptype = atypemap[this_type]
            dihs_top += [3*dih[0], 3*dih[1], 3*dih[2], this_toptype]
        else:
            dihtypes_used[this_type] = True
            ntorstypes += 1
            atypemap[this_type] = ntorstypes
            this_toptype = ntorstypes
            teq.append(RNA_params.RNA_dihs[this_type].ang)
            tk.append(RNA_params.RNA_dihs[this_type].force)
            dihs_top += [3*dih[0], 3*dih[1], 3*dih[2], this_toptype]           
    return ntorstypes,tk,teq,ndihs,dihs_top

def get_dihs(nmol, termini, CG_partnames):
    dihs = list()
    torstype = list()
    for molid in range(nmol):
        offset = termini[2*molid] - 1
        start = termini[2*molid]
        end = termini[2*molid+1] + 1
        angs_mol, torstype_mol = get_dih_list(CG_partnames[start:end])
        torstype += torstype_mol
        if offset == 0:
            dihs += angs_mol           
        else:
            for ang in angs_mol:
                at1 = ang[0]
                at2 = ang[1]
                at3 = ang[2]
                dihs.append((at1+offset, at2+offset, at3+offset))
    
    return dihs,torstype
