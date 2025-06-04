##
# @file CGmaps.py 
#
# @brief Contains the mapping information required to create the CG data from FA data
#
# Definitions to map the FA information to CG particles and then obtain the related information\n
#
# CG backbone mapping happens in the follwing way:\n 
# P -> P\n
# O5' -> O5\n
# C4' -> CA\n
# C1' -> CY\n
# O3' -> O3\n
#
# Nucleobase mapping uses the mass-weighted centre of the aromatic rigs and is as follows:\n
# A: {N7, N9, C4, C5, C8} -> A1 and {N1, C2, N3, C4, C5, N6, C6} -> A2\n
# C: {N1, C2, N3, N4, C4, C5, C6, O2} -> C1\n
# G: {N7, N9, C4, C5, C8} -> G1 and {N1, N2, C2, N3, C4, C5, C6, O6} -> G2\n
# T: {N1, C2, N3, C4, C5, C6, O2, O4} -> T1\n
# U: {N1, C2, N3, C4, C5, C6, O2, O4} -> U1\n
#
#All remaining properties for the CG particles are then defined by their name.


## Backbone map
CG_backbone_RNA = {"P": "P", "O5'": "O5", "C4'": "R4", "C1'": "R1", "O3'": "O3"}
CG_backbone_DNA = {"P": "P", "O5'": "O5", "C4'": "S4", "C1'": "S1", "O3'": "O3"}

## Nucleobase map
CG_bases = {"A1": ["N7", "N9", "C4", "C5", "C8"],
            "A2": ["N1", "C2", "N3", "C4", "C5", "N6", "C6"],
            "C1": ["N1", "C2", "N3", "N4", "C4", "C5", "C6", "O2"],
            "G1": ["N7", "N9", "C4", "C5", "C8"],
            "G2": ["N1", "N2", "C2", "N3", "C4", "C5", "C6", "O6"],
            "T1": ["N1", "C2", "N3", "C4", "C5", "C6", "O2", "O4"],
            "U1": ["N1", "C2", "N3", "C4", "C5", "C6", "O2", "O4"]}

## Defined atom masses
atom_masses = {"H": 1.008, "C": 12.011, "N": 14.007, "O": 16.000, "P": 30.974}

## CG particle masses
CG_masses = {"P": 30.970, "O5": 16.000, "C": 12.010, "R4": 20.000, "S4": 20.000, "O3": 16.000, 
             "R1": 12.000, "S1": 12.000, "A1": 67.000, "A2": 67.000, "C1": 130.000,
             "G1": 75.000, "G2": 75.000, "T1": 131.000, "U1": 131.000,
             "MG": 24.305, "NA": 22.990, "CL": 35.450}

## Reside name translation - this is redundant
CG_resnames = {"A5": "ADEi", "C5": "CYSi", "G5": "GUAi", "U5": "URAi",
               "A3": "ADE", "C3": "CYS", "G3": "GUA", "U3": "URA",
               "A": "ADE", "C": "CYS", "G": "GUA", "U": "URA",
               "DA": "ADE", "DC": "CYS", "DG": "GUA", "DT": "THY",
               "DA5": "ADEi", "DC5": "CYSi", "DG5": "GUAi", "DT5": "THYi",
               "DA3": "ADE", "DC3": "CYS", "DG3": "GUA", "DT3": "THY"}

## Simplified residue names
resnames_simple = {"A5": "A", "C5": "C", "G5": "G", "U5": "U",
                   "A3": "A", "C3": "C", "G3": "G", "U3": "U",
                   "A":  "A", "C":  "C", "G":  "G", "U":  "U",
		           "DA5": "DA", "DC5": "DC", "DG5": "DG", "DT5": "DT",
                   "DA3": "DA", "DC3": "DC", "DG3": "DG", "DT3": "DT",
                   "DA": "DA", "DC": "DC", "DG": "DG", "DT": "DT"}

# 1: O3  2: O5  3: P  4: CA  5: CY
# 6,7: G1,2  8,9: A1,2  10: U1  11: C1
# 12: D 13: MG 14: NA  15:CL

## CG particle type mapping
CG_partype = {"O3": 1, "O5": 2, "P": 3, "R4": 4,  "S4": 5, "R1": 6, "S1": 7,
              "G1": 8, "G2": 9, "A1": 10, "A2": 11, "U1": 12,
              "C1": 13, "T1": 14, "D": 15, "MG": 16, "NA": 17, "CL": 18}
# old mapping:
# CG_partype = {"C": 1, "O": 2, "P": 3, "R4": 4, "R1": 5,
#              "G1": 6, "G2": 7, "A1": 8, "A2": 9, "U1": 10,
#              "C1": 11, "D": 12, "MG": 13, "NA": 14, "CL": 15}

## CG charges mapping
CG_charges = {"O3": 0.0, "O5": 0.0, "P": -1.0, "R4": 0.0, "S4": 0.0, "R1": 0.0, "S1": 0.0, 
              "G1": 0.0, "G2": 0.0, "A1": 0.0, "A2": 0.0, "U1": 0.0,
              "C1": 0.0, "T1": 0.0, "D": 0.0, "MG": 2.0, "NA": 1.0, "CL": -1.0}
