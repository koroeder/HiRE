### CG backbone mapping ####
# P -> P
# O5' -> O5
# C5' -> C5
# C4' -> CA
# C1' -> CY

CG_backbone = {"P": "P", "O5'": "O5", "C5'": "C5", "C4'": "CA", "C1'": "CY"}

### Nucleobase mapping ####
# Mass-weighted centre of aromatic rings
# A: {N7, N9, C4, C5, C8} -> A1 and {N1, C2, N3, C4, C5, N6, C6} -> A2
# C: {N1, C2, N3, N4, C4, C5, C6, O2} -> C1
# G: {N7, N9, C4, C5, C8} -> G1 and {N1, N2, C2, N3, C4, C5, C6, O6} -> G2
# T: {N1, C2, N3, C4, C5, C6, O2, O4} -> T1
# U: {N1, C2, N3, C4, C5, C6, O2, O4} -> U1
#

CG_bases = {"A1": ["N7", "N9", "C4", "C5", "C8"],
            "A2": ["N1", "C2", "N3", "C4", "C5", "N6", "C6"],
            "C1": ["N1", "C2", "N3", "N4", "C4", "C5", "C6", "O2"],
            "G1": ["N7", "N9", "C4", "C5", "C8"],
            "G2": ["N1", "N2", "C2", "N3", "C4", "C5", "C6", "O6"],
            "T1": ["N1", "C2", "N3", "C4", "C5", "C6", "O2", "O4"],
            "U1": ["N1", "C2", "N3", "C4", "C5", "C6", "O2", "O4"]}