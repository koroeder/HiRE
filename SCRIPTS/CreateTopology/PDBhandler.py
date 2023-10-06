## @file PDBhandler.py
#
# Contains functions to deal with pdb files.\n 
# PDB file format ATOM record\n
#   Columns     Data                            Justification    Data Type\n
#    1-4         “ATOM”                          character\n
#    7-11        Atom serial number              right            integer\n
#    13-16       Atom name                       left             character\n
#    17          Alternate location indicator                     character\n
#    18-20       Residue name                    right            character\n
#    22          Chain identifier                                 character\n
#    23-26       Residue sequence number         right            integer\n
#    27          Code for residue insertion                       character\n
#    31-38       X orthogonal Å coordinate       right            real (8.3)\n
#    39-46       Y orthogonal Å coordinate       right            real (8.3)\n
#    47-54       Z orthogonal Å coordinate       right            real (8.3)\n
#    55-60       Occupancy                       right            real (6.2)\n
#    61-66       Temperature factor              right            real (6.2)\n
#    73-76       Segment identifier              left             character\n
#    77-78       Element symbol                  right            character\n

import sys

## Function to parse line in pdb
#
# The parsing is based on the right alignment as detailed in the file descriptions. 
# The element might be empty depending on the file format.
#
# @return Returns the atom name, residue name, residue id, atom coordinates and element for each entry
def parse_line(line):
    atom_name = line[12:16].strip()
    res_name = line[17:20].strip()
    res_id = int(line[22:26])
    atom_xyz = [float(line[30:38]), float(line[38:46]), float(line[46:54])]
    element = line[76:78].strip()
    return [atom_name, res_name, res_id, atom_xyz, element]

## Function to open pdb and parse it line by line
#
# Returns the atom data as dictionary,
# the number of atoms (natom) and a list of terminal atoms 
# Recognition of termini is based on TER statements in pdb file
#
# @param[in] inpfile - input file name
#
# @return Dictionary containg all pdb data, numebr of atoms and a list of terminal atoms
def get_pdb_data(inpfile):
    usemodels = False
    data = dict()
    termini = [1] #first atom is always a terminus
    natom = 0
    with open(inpfile,"r") as f:
        lines = f.readlines()
        for line in lines:
            if line[:5] == "MODEL":
                if usemodels:
                    return data,termini,natom
                else:
                    usemodels = True
            #parse atom information
            if line[:4] == "ATOM":
                natom += 1
                data[natom] = parse_line(line)
            #add entry for every terminal atom
            elif line[:3] == "TER":
                termini.append(natom)
            else:
                continue
    return data,termini,natom

## Obtain residue information
#
# Function to obtain the number of residues (resid, dual used as counter),
# the first and last atom of each residue as dictionary (res), a list of
# all residue names (resnames), and the coordinates for each residue (coordsbyres).
#
# @param[in] natom - Number of atoms
# @param[in] data - pdb information from parsing
def get_residues(natom,data):
    res = {1: [1,0]}
    resid = 1
    for idx in range(natom):
        if resid < data[idx+1][2]:
            res[resid][1] = idx
            resid += 1
            res[resid] = [idx+1,0]
    res[resid][1] = natom
    resnames = list()
    for idx in range(resid):
        resnames.append(data[res[idx+1][0]][1])
    coordsbyres = dict()
    for ridx in range(resid):
        coords = list()
        for idx in range(res[ridx+1][0],res[ridx+1][1]+1):
            coords.append(data[idx][3])
        coordsbyres[ridx+1] = coords
    return resid, res, resnames, coordsbyres

def shift_resids(natom,data):
    resid = data[1][2]
    shift = resid - 1
    prev = 1
    for i in range(natom):
        curr = data[i+1][2] - shift
        if prev<curr:
            if prev+1==curr:
                prev = curr
            else:
                shift = shift + curr - prev - 1
                curr = data[i+1][2] - shift
                prev = curr
        data[i+1][2] = curr
    return data

## Function to get a list of elements from data dictionary
def get_elements(natom,data):
    elements = list()
    for idx in range(natom):
        elements.append(data[idx+1][4])
    return elements

## Function to get a list of atom names from data dictionary
def get_atomnames(natom,data):
    atomnames = list()
    for idx in range(natom):
        atomnames.append(data[idx+1][0])
    return atomnames

## Function to make sure we have correct termini
#
# The final entry should be the last atom
# In addition, we added the last atom of each molecule, but not the first
# Here, we fix this problem, and then transfer the list into a set to
# remove duplicates and then back into a sorted list to be useful
def fix_termini(natom,termini):
    if termini[-1] != natom:
        termini.append(natom)
    new_term = list()
    for t in termini:
        if t==1 or t==natom:
            continue
        else:
            new_term.append(t+1)
    return sorted(set(termini+new_term))

## Function to make sure we have correct termini
#
# The final entry should be the last atom
# In addition, we added the last atom of each molecule, but not the first
# Here, we fix this problem, and then transfer the list into a set to
# remove duplicates and then back into a sorted list to be useful
def fix_termini(natom,termini):
    if termini[-1] != natom:
        termini.append(natom)
    new_term = list()
    for t in termini:
        if t==1 or t==natom:
            continue
        else:
            new_term.append(t+1)
    return sorted(set(termini+new_term))

def fix_termini_CG(natom,termini):
    if termini[-1] != natom:
        termini.append(natom)
    new_term = list()
    for t in termini:
        if t==1 or t==natom:
            continue
        else:
            new_term.append(t+1)
    return sorted(termini+new_term)

## Function to parse file and obtain all relevant data
def parse_pdb(inpfile):
    data, termini, natom = get_pdb_data(inpfile)
    data = shift_resids(natom,data)
    nres, res, resnames, coordsbyres= get_residues(natom,data)
    elements = get_elements(natom,data)
    atomnames = get_atomnames(natom,data)
    if elements[0] == "":
        for idx,name in enumerate(atomnames):
            elements[idx] = name[0]
    termini = fix_termini(natom, termini)
    return natom,nres,atomnames,elements,res,resnames,coordsbyres,termini

## Function to parse CG file and obtain all relevant data
def parse_CGpdb(inpfile):
    data, termini, natom = get_pdb_data(inpfile)
    data = shift_resids(natom,data)
    nres, res, resnames, coordsbyres= get_residues(natom,data)
    elements = get_elements(natom,data)
    atomnames = get_atomnames(natom,data)
    if elements[0] == "":
        for idx,name in enumerate(atomnames):
            elements[idx] = name[0]
    termini = fix_termini_CG(natom, termini)
    return natom,nres,atomnames,elements,res,resnames,coordsbyres,termini

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


# Can run this script alone to check everything is parsed correctly
if __name__ == "__main__":
    natom,nres,atomnames,elements,res,resnames,coordsbyres,termini = parse_pdb(sys.argv[1])
    print("Number of atoms: ", natom)
    print("Number of residues: ", nres)
    print("Atom names: ", atomnames)
    print("Elements: ", elements)
    print("Residues: ", res)
    print("Residue names: ", resnames)
    print("Termini: ", termini)