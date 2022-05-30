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

