import CGmaps

def massW_centre(coords,masses):
    X = 0.0
    Y = 0.0
    Z = 0.0
    mtot = 0.0
    for idx,m in enumerate(masses):
        X += m*coords[idx][0]
        Y += m*coords[idx][1]
        Z += m*coords[idx][2]
        mtot += m
    return [X/mtot,Y/mtot,Z/mtot]

def assign_atommass(atomnames):
    atommasses = list()
    for name in atomnames:
        atommasses.append(CGmaps.atom_masses[name[0]])
    return atommasses

def assign_CGmasses(CGnames):
    CGmasses = list()
    for name in CGnames:
        CGmasses.append(CGmaps.CG_masses[name])
    return CGmasses

def FA2CG_res(resid,labels,masses,coords):
    if (resid == "A" or resid == "DA"):
        Part_labels = CGmaps.CG_bases["A1"] + CGmaps.CG_bases["A2"]
        name_p1 = "A1"  
        name_p2 = "A2"        
        parts = 2
    elif (resid == "C" or resid == "DC"):
        Part_labels = CGmaps.CG_bases["C1"]
        name_p1 = "C1"        
        parts = 1
    elif (resid == "G" or resid == "DG"):
        Part_labels = CGmaps.CG_bases["G1"] + CGmaps.CG_bases["G2"]
        name_p1 = "G1"  
        name_p2 = "G2"                  
        parts = 2
    elif (resid == "DT"):
        Part_labels = CGmaps.CG_bases["T1"]
        name_p1 = "T1"       
        parts = 1
    elif (resid == "U"):
        Part_labels = CGmaps.CG_bases["U1"]
        name_p1 = "U1"
        parts = 1

    CG_labels = list()
    CG_coords = list()
    part_mass = list()
    part_coords = list()
    part_mass2 = list()
    part_coords2 = list()

    for idx,l in enumerate(labels):
        if l in CGmaps.CG_backbone.keys():
            CG_labels.append(CGmaps.CG_backbone[l])
            CG_coords.append(coords[idx][:])
        elif l in Part_labels:
            if parts == 1:
                part_mass.append(masses[idx])
                part_coords.append(coords[idx][:])
            else:
                if ((l in CGmaps.CG_bases["A1"]) or (l in CGmaps.CG_bases["G1"])):
                    part_mass.append(masses[idx])
                    part_coords.append(coords[idx][:])
                if ((l in CGmaps.CG_bases["A2"]) or (l in CGmaps.CG_bases["G2"])):
                    part_mass2.append(masses[idx])
                    part_coords2.append(coords[idx][:])             
    if parts == 1:
        CG_labels.append(name_p1)
        CG_coords.append(massW_centre(part_coords,part_mass))
    else:
        CG_labels.append(name_p1)
        CG_coords.append(massW_centre(part_coords,part_mass))       
        CG_labels.append(name_p2)
        CG_coords.append(massW_centre(part_coords2,part_mass2))
    
    return CG_labels, CG_coords


def get_CG_info(reslabel,resrange,atomlabels,masses,coords):
    CG_labels = list()
    CG_coords = list()
    CG_resrange = dict()
    start_p = 1
    end_p = 0
    for idx,label in enumerate(reslabel):
        start = resrange[idx+1][0]
        end = resrange[idx+1][1]
        res_masses = masses[start:end+1]
        res_coords = coords[idx+1]
        res_labels = atomlabels[start:end+1]
        resCG_labels, resCG_coords = FA2CG_res(label, res_labels, res_masses, res_coords)
        CG_labels += resCG_labels
        CG_coords += resCG_coords
        start_p = end_p + 1
        end_p = end_p + len(resCG_labels)
        CG_resrange[idx+1] = [start_p,end_p]
    return CG_labels,CG_coords,CG_resrange

def get_CG_termini(resrange, termini, CG_resrange):
    CG_termini = list()
    CG_teri_flags = list()
    for idx,ids in enumerate(resrange):
        start = ids[0]
        end = ids[1]
        if start in termini:
            CG_termini.append(CG_resrange[idx][0])
            CG_teri_flags.append(True)
        if end in termini:
            CG_termini.append(CG_resrange[idx][1])
        if len(CG_teri_flags)<(idx+1):
            CG_teri_flags.append(False)
    return CG_termini,CG_teri_flags           

def assign_resnames(resnames):
    CG_resnames = list()
    new_resnames = list()
    for name in resnames:
        CG_resnames.append(CGmaps.CG_resnames[name])
        new_resnames.append(CGmaps.resnames_simple[name])
    return new_resnames, CG_resnames

def assign_par_properties(CG_labels):
    CG_charges = list()
    CG_partype = list()
    for name in CG_labels:
        CG_charges.append(CGmaps.CG_charges[name])
        CG_partype.append(CGmaps.CG_partype[name])
    return CG_charges, CG_partype