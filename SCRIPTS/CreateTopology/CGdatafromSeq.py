## @file CGdatafromSeq.py
#
# Contains functions to get cG information from sequence
import CGmaps


def find_termini_seqlist(seqlist):
## Find termini residue ids from sequence
#
# param[in] seqlist - List of residue names 
    CGres_termini = list()
    for idx,res in enumerate(seqlist):
        if res in ["A5", "C5", "G5", "U5", "DA5", "DC5", "DG5", "DT5"]:
            CGres_termini.append(idx+1)
        elif res in ["A3", "C3", "G3", "U3", "DA3", "DC3", "DG3", "DT3"] :
            CGres_termini.append(idx+1)
        else:
            continue
    return CGres_termini

def get_atoms_from_seq(seqlist):
## Get list of CG particle names from sequence
#
# param[in] seqlist - List of residue names
    CG_atoms = list()
    for res in seqlist:
        if res in ["A5", "C5", "G5", "U5", "DA5", "DC5", "DG5", "DT5"]:
            bb = ["O", "C", "R4", "R1"]
        else:
            bb = ["P", "O", "C", "R4", "R1"]
        if res in ["A3", "A5", "A", "DA", "DA3", "DA5"]:
            sc = ["A1", "A2"]
        elif res in ["C3", "C5", "C", "DC", "DC3", "DC5"]:
            sc = ["C1"]
        elif res in ["G3", "G5", "G", "DG", "DG3", "DG5]:
            sc = ["G1", "G2"]
        elif res in ["DT3", "DT5", "DT"]:
            sc = ["T1"]
        elif res in ["U3", "U5", "U"]:
            sc = ["U1"]                  
        CG_atoms = CG_atoms + bb + sc
    return CG_atoms

 
def get_resrange(seqlist):
## Get list of id ranges for each residue
#
# param[in] seqlist - List of residue names
    nstart = 1
    CG_resrange = dict()
    for res in seqlist:
        if res in ["C5", "U5", "DC5", "DT5"]:
            size = 5
        elif res in ["A5", "G5", "DA5", "DG5"]:
            size = 6
        elif res in ["C", "C3", "U", "U3", "DC", "DC3", "DT3", "DT"]:
            size = 6
        elif res in ["A", "A3", "G", "G3", "DA", "DA3", "DG", "DG3"]:
            size = 7
        CG_resrange[len(CG_resrange)+1] = [nstart, nstart + size -1]
        nstart = nstart + size
    return CG_resrange

def get_CGterminis(CG_resrange,termini):
## Get CG termini particle indices
#
# param[in] CG_resrange - dictionary of first and last particle in each residue
# param[in] termini - list of terminal residues (3' and 5' nucleotides)
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


def translate_CGresnames(seqlist):
## Return list of new residue names
# TODO: remove this function as it is not needed
    new_resnames = list()
    for res in seqlist:
        new_resnames.append(CGmaps.CG_resnames[res])
    return new_resnames
