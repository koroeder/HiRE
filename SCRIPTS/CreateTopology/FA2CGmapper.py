## @file FA2CGmapper.py
#
# Set of functions to obtain HiRE coarse-grained information from atom data

import CGmaps

## Function to get centre of mass coordinates
#
# @param[in] coords - coordinates of particles
# @param[in] masses - masses of particles
#
# @return List with centre of mass coordinates
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

## Function to assign atom masses
#
# As the first letter of the atom names in pdbs generally correspond to
# the element, and the element information in the pdb might be empty,
# we use the first letter of the atom name to map to the atom mass.
# The map is defined in CGmaps.py.
#
# @param[in] atomnames - names of atoms
#
# @return List of atom masses
def assign_atommass(atomnames):
    atommasses = list()
    for name in atomnames:
        atommasses.append(CGmaps.atom_masses[name[0]])
    return atommasses

## Function to assign the CG particle mass based on particle name
# The map is defined in CGmaps.py.
#
# @param[in] CGnames - names of CG particles
#
# @return List of CG particle masses
def assign_CGmasses(CGnames):
    CGmasses = list()
    for name in CGnames:
        CGmasses.append(CGmaps.CG_masses[name])
    return CGmasses

## Function to map FA to CG for a single residue
#
# @param[in] resid - index of residue
# @param[in] labels - list of atom names
# @param[in] masses - list of atom masses
# @param[in] coords - coordinates
#
# @return A list of CG particle names and the corresponding CG coordinates
def FA2CG_res(resid,labels,masses,coords):
    # First we define the atoms we are interest in
    # This will depend on the nucleotide
    # Importantly, we obtain the dummy variables for the number of base particles
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

    # Iterate over every atom in the residue, and follow the appropriate CG mapping
    for idx,l in enumerate(labels):
        # For all heavy backbone atoms, we map 1:1 as given in CGmaps
        if l in CGmaps.CG_backbone.keys():
            CG_labels.append(CGmaps.CG_backbone[l])
            CG_coords.append(coords[idx][:])
        # For any atom in the base particles, we append the relevant coordinates 
        # and mass to a list. Here, we distinguish between a 
        # single particle and two particles
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
    # Now we obtain the coordinates and mass for the base particles
    if parts == 1:
        CG_labels.append(name_p1)
        CG_coords.append(massW_centre(part_coords,part_mass))
    else:
        CG_labels.append(name_p1)
        CG_coords.append(massW_centre(part_coords,part_mass))       
        CG_labels.append(name_p2)
        CG_coords.append(massW_centre(part_coords2,part_mass2))
    
    return CG_labels, CG_coords

## Function to get CG information
# The three sets of data are the particle names, their coordinates 
# and the first and last particle for each residue
#
# @param[in] reslabel - List of residue names
# @param[in] resrange - Dictionary containing the first and last atom index for each residue
# @param[in] atomlabels - List of atom names
# @param[in] masses - List of atom masses
# @param[in] coords - Coordinates
#
# @return List of CG particle labels, CG coordinates and dictionary of first and last particle for each residue
def get_CG_info(reslabel,resrange,atomlabels,masses,coords):
    CG_labels = list()
    CG_coords = list()
    CG_resrange = dict()
    start_p = 1
    end_p = 0
    for idx,label in enumerate(reslabel):
        start = resrange[idx+1][0]
        end = resrange[idx+1][1]
        res_masses = masses[start-1:end]
        res_coords = coords[idx+1]
        res_labels = atomlabels[start-1:end]
        resCG_labels, resCG_coords = FA2CG_res(label, res_labels, res_masses, res_coords)
        CG_labels += resCG_labels
        CG_coords += resCG_coords
        start_p = end_p + 1
        end_p = end_p + len(resCG_labels)
        CG_resrange[idx+1] = [start_p,end_p]
    return CG_labels,CG_coords,CG_resrange

## @brief Obtain CG termini
#
# Function to produce the terminal atoms for the molecules in our system.\n
# This information is obtained by mapping the termini for the atomistic 
# system to the CG system using the dictionaries with the first and last 
# atom/particle for each residue.
#
# @param[in] resrange - Dictionary of residue ranges
# @param[in] termini - Information on terminal residues
# @param[in] resrange - Dictionary of residue ranges for CG
#
# @return CG terminal information
def get_CG_termini(resrange, termini, CG_resrange):
    CG_termini = list()
    for key in resrange.keys():
        start = resrange[key][0]
        end = resrange[key][1]
        if start in termini:
            CG_termini.append(CG_resrange[key][0])
        if end in termini:
            CG_termini.append(CG_resrange[key][1])
    return CG_termini         

## Redundant apart form simplified information!
#
# Here, we assign the proper CG residue names, including terminal information
# We also produce a simplified version of the residue names,
# which is used in the functions that assign bonds, angles etc.
# The maps are again defined in CGmaps.py
def assign_resnames(resnames):
    CG_resnames = list()
    new_resnames = list()
    for name in resnames:
        CG_resnames.append(CGmaps.CG_resnames[name])
        new_resnames.append(CGmaps.resnames_simple[name])
    return new_resnames, CG_resnames

## Function to assign CG charges and particle type
# Again, we use fixed maps for this purpose and use the particle name
def assign_par_properties(CG_labels):
    CG_charges = list()
    CG_partype = list()
    for name in CG_labels:
        CG_charges.append(CGmaps.CG_charges[name])
        CG_partype.append(CGmaps.CG_partype[name])
    return CG_charges, CG_partype
