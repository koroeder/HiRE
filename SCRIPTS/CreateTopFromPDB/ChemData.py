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



