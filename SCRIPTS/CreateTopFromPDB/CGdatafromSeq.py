import CGmaps

# find termini residue ids from sequence 
def find_termini_seqlist(seqlist):
    CGres_termini = list()
    for idx,res in enumerate(seqlist):
        if res in ["A5", "C5", "G5", "T5", "U5"]:
            CGres_termini.append(idx+1)
        elif res in ["A3", "C3", "G3", "T3", "U3"] :
            CGres_termini.append(idx+1)
        else:
            continue
    return CGres_termini

# get list of atoms from sequence
def get_atoms_from_seq(seqlist):
    CG_atoms = list()
    for res in seqlist:
        if res in ["A5", "C5", "G5", "T5", "U5"]:
            bb = ["O", "C", "R4", "R1"]
        else:
            bb = ["P", "O", "C", "R4", "R1"]
        if res in ["A3", "A5", "A", "DA"]:
            sc = ["A1", "A2"]
        elif res in ["C3", "C5", "C", "DC"]:
            sc = ["C1"]
        elif res in ["G3", "G5", "G", "DG"]:
            sc = ["G1", "G2"]
        elif res in ["T3", "T5", "DT"]:
            sc = ["T1"]
        elif res in ["U3", "U5", "U"]:
            sc = ["U1"]                  
        # now add all entries to the list
        CG_atoms = CG_atoms + bb + sc
    return CG_atoms

# get list of ranges for each residue
def get_resrange(seqlist):
    nstart = 1
    CG_resrange = dict()
    for res in seqlist:
        if res in ["C5", "T5", "U5"]:
            size = 5
        elif res in ["A5", "G5"]:
            size = 6
        elif res in ["C3", "T3", "U3", "C", "U", "DC", "DT"]:
            size = 6
        elif res in ["A3", "G3", "A", "G", "DA", "DG"]:
            size = 7
        CG_resrange[len(CG_resrange)+1] = [nstart, nstart + size -1]
        nstart = nstart + size
    return CG_resrange

#get CG termini particle indices
def get_CGterminis(CG_resrange,termini):
    CG_termini = list()
    start = True
    for res in CG_resrange.keys():
        if res in termini:
            if start:
                CG_termini.append(CG_resrange[res][0])
                start = False
            else:
                CG_termini.append(CG_resrange[res][1])
                start = True
    return CG_termini              


#return list of new residue names
def translate_CGresnames(seqlist):
    new_resnames = list()
    for res in seqlist:
        new_resnames.append(CGmaps.CG_resnames[res])
    return new_resnames
