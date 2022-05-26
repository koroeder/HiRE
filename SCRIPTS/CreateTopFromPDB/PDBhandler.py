#### PDB file format ATOM record
#	Columns	Data			Justification	Data Type
#	1-4	“ATOM”						character
#	7-11	Atom serial number		right		integer
#	13-16	Atom name			left		character
#	17	Alternate location indicator			character
#	18-20	Residue name			right		character
#	22	Chain identifier				character
#	23-26	Residue sequence number	right		integer
#	27	Code for residue insertion			character
#	31-38	X orthogonal Å coordinate	right		real (8.3)
#	39-46	Y orthogonal Å coordinate	right		real (8.3)
#	47-54	Z orthogonal Å coordinate	right		real (8.3)
#	55-60	Occupancy			right		real (6.2)
#	61-66	Temperature factor		right		real (6.2)
#	73-76	Segment identifier		left		character
#	77-78	Element symbol			right		character

import sys


def parse_line(line):
    atom_name = line[12:16].strip()
    res_name = line[17:20].strip()
    res_id = int(line[22:26])
    atom_xyz = [float(line[30:38]), float(line[38:46]), float(line[46:54])]
    element = line[76:78].strip()
    return atom_name, res_name, res_id, atom_xyz, element

def get_pdb_data(inpfile):
    data = dict()
    termini = [1] #first atom is always a terminus
    natom = 0
    with open(inpfile,"r") as f:
        lines = f.readlines()
        for line in lines:
            if line[:4] == "ATOM":
                natom += 1
                data[natom] = parse_line(line)
            elif line[:3] == "TER":
                termini.append(natom)
            else:
                continue
    return data,termini,natom

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

def get_elements(natom,data):
    elements = list()
    for idx in range(natom):
        elements.append(data[idx+1][4])
    return elements

def get_atomnames(natom,data):
    atomnames = list()
    for idx in range(natom):
        atomnames.append(data[idx+1][0])
    return atomnames

def parse_pdb(inpfile):
    data, termini, natom = get_pdb_data(inpfile)
    nres, res, resnames, coordsbyres= get_residues(natom,data)
    elements = get_elements(natom,data)
    atomnames = get_atomnames(natom,data)
    if elements[0] == "":
        for idx,name in enumerate(atomnames):
            elements[idx] = name[0]
    return natom,nres,atomnames,elements,res,resnames,coordsbyres,termini

if __name__ == "__main__":
    natom,nres,atomnames,elements,res,resnames,coordsbyres,termini = parse_pdb(sys.argv[1])
    print("Number of atoms: ", natom)
    print("Number of residues: ", nres)
    print("Atom names: ", atomnames)
    print("Elements: ", elements)
    print("Residues: ", res)
    print("Residue names: ", resnames)
    print("Termini: ", termini)