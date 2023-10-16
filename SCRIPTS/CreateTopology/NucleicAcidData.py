#Classes and data to define Nucleic acids
#Currently not used
class Atom_NA:
    def __init__(self, idx, name, el, mass):
        self.idx = idx
        self.name = name
        self.el = el
        self.mass = mass

    def set_idx(self,newidx):
        self.idx = newidx

class Bond_NA:
    def __init__(self, typeid, atom1, atom2, force, dist):
        self.typeid = typeid
        self.atom1 = atom1
        self.atom2 = atom2
        self.force = force
        self.dist = dist
    
    def iscorrectbond(self, at1, at2):
        if ((at1 is self.atom1) and (at2 is self.atom2)):
            return True
        elif ((at2 is self.atom1) and (at1 is self.atom2)):
            return True
        else:
            return False
    
    def get_typeid(self):
        return self.typeid

class Angle_NA:
    def __init__(self, name, typeid, atom1, atom2, atom3, force, ang):
        self.name = name
        self.typeid = typeid
        self.atom1 = atom1
        self.atom2 = atom2
        self.atom3 = atom3
        self.force = force
        self.ang = ang
        
    def iscorrectangle(self, at1, at2, at3):
        if ((at1 is self.atom1) and (at2 is self.atom2) and (at3 is self.atom3)):
            return True
        else:
            return False
        
    def get_typeid(self):
        return self.typeid
        
class Dihedral_NA:
    def __init__(self, name, typeid, atom1, atom2, atom3, atom4, force, dih, nterm):
        self.name = name
        self.typeid = typeid
        self.atom1 = atom1
        self.atom2 = atom2
        self.atom3 = atom3
        self.atom4 = atom4
        self.force = force
        self.dih = dih
        self.nterm = nterm
    
    def iscorrectdih(self, at1, at2, at3, at4):
        if ((at1 is self.atom1) and (at2 is self.atom2) and (at3 is self.atom3) and (at4 is self.atom4)):
            return True
        else:
            return False  

    def get_typeid(self):
        return self.typeid	      

RA = [	Atom_NA(1,   "P",  "P",   36.970),
	Atom_NA(2, "O5*",  "O",   16.000),
	Atom_NA(3, "C5*",  "C",   12.010),
	Atom_NA(4,  "CA", "R4",   20.000),
	Atom_NA(5,  "CY", "R1",   12.000),
	Atom_NA(6,  "A1", "A1",   67.000),
	Atom_NA(7,  "A2", "A2",   67.000)]

DA = [	Atom_NA(1,   "P",  "P",   36.970),
	Atom_NA(2, "O5*",  "O",   16.000),
	Atom_NA(3, "C5*",  "C",   12.010),
	Atom_NA(4,  "CA", "R4",   20.000),
	Atom_NA(5,  "CY", "S1",   12.000),
	Atom_NA(6,  "A1", "A1",   67.000),
	Atom_NA(7,  "A2", "A2",   67.000)]

RC = [	Atom_NA(1,   "P",  "P",   36.970),
	Atom_NA(2, "O5*",  "O",   16.000),
	Atom_NA(3, "C5*",  "C",   12.010),
	Atom_NA(4,  "CA", "R4",   20.000),
	Atom_NA(5,  "CY", "R1",   12.000),
	Atom_NA(6,  "C1", "C1",  130.000)]

DC = [	Atom_NA(1,   "P",  "P",   36.970),
	Atom_NA(2, "O5*",  "O",   16.000),
	Atom_NA(3, "C5*",  "C",   12.010),
	Atom_NA(4,  "CA", "R4",   20.000),
	Atom_NA(5,  "CY", "S1",   12.000),
	Atom_NA(6,  "C1", "C1",  130.000)]

RG = [	Atom_NA(1,   "P",  "P",   36.970),
	Atom_NA(2, "O5*",  "O",   16.000),
	Atom_NA(3, "C5*",  "C",   12.010),
	Atom_NA(4,  "CA", "R4",   20.000),
	Atom_NA(5,  "CY", "R1",   12.000),
	Atom_NA(6,  "G1", "G1",   75.000),
	Atom_NA(7,  "G2", "G2",   75.000)]

DG = [	Atom_NA(1,   "P",  "P",   36.970),
	Atom_NA(2, "O5*",  "O",   16.000),
	Atom_NA(3, "C5*",  "C",   12.010),
	Atom_NA(4,  "CA", "R4",   20.000),
	Atom_NA(5,  "CY", "S1",   12.000),
	Atom_NA(6,  "G1", "G1",   75.000),
	Atom_NA(7,  "G2", "G2",   75.000)]

DT = [	Atom_NA(1,   "P",  "P",   36.970),
	Atom_NA(2, "O5*",  "O",   16.000),
	Atom_NA(3, "C5*",  "C",   12.010),
	Atom_NA(4,  "CA", "R4",   20.000),
	Atom_NA(5,  "CY", "S1",   12.000),
	Atom_NA(6,  "T1", "T1",  131.000)]

RU = [	Atom_NA(1,   "P",  "P",   36.970),
	Atom_NA(2, "O5*",  "O",   16.000),
	Atom_NA(3, "C5*",  "C",   12.010),
	Atom_NA(4,  "CA", "R4",   20.000),
	Atom_NA(5,  "CY", "R1",   12.000),
	Atom_NA(6,  "U1", "U1",  131.000)]

NA_atoms = [	Atom_NA(0,   "P",  "P",   36.970),
		Atom_NA(0, "O5*",  "O",   16.000),
		Atom_NA(0, "C5*",  "C",   12.010),
		Atom_NA(0,  "CA", "R4",   20.000),
		Atom_NA(0,  "CY", "R1",   12.000),  # RNA sugar
 		Atom_NA(0,  "CY", "S1",   12.000),  # DNA sugar     
		Atom_NA(0,  "A1", "A1",   67.000),
		Atom_NA(0,  "A2", "A2",   67.000),
		Atom_NA(0,  "G1", "G1",   75.000),
		Atom_NA(0,  "G2", "G2",   75.000),
		Atom_NA(0,  "C1", "C1",  130.000),
        Atom_NA(0,  "T1", "T1",  131.000),
		Atom_NA(0,  "U1", "U1",  131.000),
		Atom_NA(0,   "D",  "D", 1000.000)]

NA_bonds = [	
        Bond_NA( 1, "R4",  "P",  30, 3.800),
		Bond_NA( 2, "R4", "R1", 200, 2.344), #RNA
		Bond_NA( 3, "R4", "S1", 200, 2.344), #DNA   
		Bond_NA( 4, "R1", "G1", 200, 2.622), #RNA
		Bond_NA( 5, "R1", "A1", 200, 2.633), #RNA
		Bond_NA( 6, "R1", "U1", 200, 3.062), #RNA
		Bond_NA( 7, "R1", "C1", 200, 3.004), #RNA
		Bond_NA( 8, "S1", "G1", 200, 2.622), #DNA
		Bond_NA( 9, "S1", "A1", 200, 2.633), #DNA
		Bond_NA(10, "S1", "T1", 200, 3.062), #DNA
		Bond_NA(11, "S1", "C1", 200, 3.004), #DNA   
		Bond_NA(12, "G1", "G2", 200, 2.450),
		Bond_NA(13, "A1", "A2", 200, 2.180),
		Bond_NA(14,  "C", "R4", 200, 1.520),
		Bond_NA(15,  "P",  "O", 200, 1.593),
		Bond_NA(16,  "O",  "C", 200, 1.430),
		Bond_NA(17,  "D",  "D",  40, 12.00),
		Bond_NA(18, "R4",  "D",  10, 12.00),
		Bond_NA(19,  "D",  "O",  10, 12.00),
		Bond_NA(20,  "D",  "P",  10, 12.00)]

# For DNA: Bond_DNA("R4",  "P",  30, 3.950

NA_angles = [	
        Angle_NA("RNA01",  1, "R4", "R1", "A1",  70, 123.6),
		Angle_NA("RNA02",  2, "R4", "R1", "U1",  70, 132.5),
		Angle_NA("RNA03",  3, "R4", "R1", "G1",  70, 123.5),
		Angle_NA("RNA04",  4, "R4", "R1", "C1",  70, 131.1),
		Angle_NA("RNA05",  5, "R1", "A1", "A2", 120, 116.7),
		Angle_NA("RNA06",  6, "R1", "G1", "G2", 120, 111.2),
		Angle_NA("RNA07",  7,  "P",  "O",  "C",  70, 122.9),
		Angle_NA("RNA08",  8,  "O",  "C", "R4",  70, 110.6),
		Angle_NA("RNA09",  9,  "C", "R4",  "P",  70,  98.0),
		Angle_NA("RNA10", 10, "R4",  "P",  "O",  50, 110.0),    
		Angle_NA("RNA11", 11,  "C", "R4", "R1",  70, 135.5),
		Angle_NA("RNA12", 12, "R1", "R4",  "P", 100,  98.0),
		Angle_NA("RNA13", 13,  "D",  "D",  "D",  80, 180.0),
        Angle_NA("DNA01", 14, "R4", "S1", "A1",  70, 130.0),
		Angle_NA("DNA02", 15, "R4", "S1", "T1",  70, 142.0),
		Angle_NA("DNA03", 16, "R4", "S1", "G1",  70, 131.0),
		Angle_NA("DNA04", 17, "R4", "S1", "C1",  70, 146.0),
		Angle_NA("DNA05", 18, "S1", "A1", "A2", 120, 116.0),
		Angle_NA("DNA06", 19, "S1", "G1", "G2", 120, 111.2),
		Angle_NA("DNA07", 20,  "P",  "O",  "C",  70, 120.0),
		Angle_NA("DNA08", 21,  "O",  "C", "R4",  70, 110.6),
		Angle_NA("DNA09", 22,  "C", "R4",  "P",  70, 110.0),
		Angle_NA("DNA10", 23, "R4",  "P",  "O",  50, 105.0),    
		Angle_NA("DNA11", 24,  "C", "R4", "S1",  70, 130.0),
		Angle_NA("DNA12", 25, "S1", "R4",  "P", 100,  85.0),
		Angle_NA("DNA13", 26,  "D",  "D",  "D",  80, 180.0)]

NA_dihs = [	
        Dihedral_NA("RNA01",  1, "R4", "R1", "G1", "G2", 1.0,     -20, 1),
		Dihedral_NA("RNA02",  2, "R4", "R1", "G1", "G2", 0.2,    -160, 1),
		Dihedral_NA("RNA03",  3, "R4", "R1", "A1", "A2", 1.0,     -20, 1),
		Dihedral_NA("RNA04",  4, "R4", "R1", "A1", "A2", 0.2,  -160.0, 1),    
		Dihedral_NA("RNA05",  5, "R4", "A1", "A2", "R1", 1.0,  -165.0, 1),
		Dihedral_NA("RNA06",  6, "R4", "G1", "G2", "R1", 1.0,  -165.0, 1),
		Dihedral_NA("RNA07",  7,  "C", "R4", "R1", "A1", 1.0,  -150.0, 1),
		Dihedral_NA("RNA08",  8,  "C", "R4", "R1", "A1", 0.2,   150.0, 1),
		Dihedral_NA("RNA09",  9,  "C", "R4", "R1", "G1", 1.0,  -150.0, 1),
		Dihedral_NA("RNA10", 10,  "C", "R4", "R1", "C1", 1.0,  -150.0, 1),
		Dihedral_NA("RNA11", 11,  "C", "R4", "R1", "C1", 0.2,   150.0, 1),
		Dihedral_NA("RNA12", 12,  "C", "R4", "R1", "U1", 1.0,  -150.0, 1),
		Dihedral_NA("RNA13", 13,  "C", "R4", "R1", "U1", 0.2,   150.0, 1),
		Dihedral_NA("RNA14", 14,  "P", "R4", "R1", "A1", 1.0,   100.0, 1),
		Dihedral_NA("RNA15", 15,  "P", "R4", "R1", "G1", 1.0,   100.0, 1),
		Dihedral_NA("RNA16", 16,  "P", "R4", "R1", "C1", 1.0,   100.0, 1),
		Dihedral_NA("RNA17", 17,  "P", "R4", "R1", "U1", 1.0,   100.0, 1),
		Dihedral_NA("RNA18", 18,  "C", "R4",  "P",  "O", 1.2,    20.0, 1),
		Dihedral_NA("RNA19", 19, "R1", "R4",  "P",  "O", 1.2,   160.0, 1),
		Dihedral_NA("RNA20", 20,  "O",  "C", "R4",  "P", 1.0,    90.0, 3),
		Dihedral_NA("RNA21", 21,  "O",  "C", "R4", "R1", 1.0,   150.0, 1),
		Dihedral_NA("RNA22", 22,  "P",  "O",  "C", "R4", 0.33,    0.0, 1),
		Dihedral_NA("RNA23", 23,  "P",  "O",  "C", "R4", 0.125, 180.0, 2),
		Dihedral_NA("RNA24", 24,  "P",  "O",  "C", "R4", 0.83,    0.0, 3),
		Dihedral_NA("RNA25", 25, "R4",  "P",  "O",  "C", 1.00,    0.0, 2),
        
        Dihedral_NA("DNA01", 26, "R4", "S1", "G1", "G2", 1.0,  -120.0, 2),
		Dihedral_NA("DNA02", 27, "R4", "S1", "A1", "A2", 1.0,  -120.0, 2),
		Dihedral_NA("DNA03", 28, "R4", "S1", "A1", "A2", 0.2,   160.0, 1),    
		Dihedral_NA("DNA04", 29, "R4", "A1", "A2", "S1", 1.0,   165.0, 1),
		Dihedral_NA("DNA05", 30,  "C", "R4", "S1", "A1", 1.0,   150.0, 1),
		Dihedral_NA("DNA06", 31,  "C", "R4", "S1", "G1", 1.0,   150.0, 1),
		Dihedral_NA("DNA07", 32,  "C", "R4", "S1", "C1", 1.0,   150.0, 1),
		Dihedral_NA("DNA08", 33,  "C", "R4", "S1", "U1", 1.0,   150.0, 1),
		Dihedral_NA("DNA09", 34,  "P", "R4", "S1", "A1", 1.0,    30.0, 1),
		Dihedral_NA("DNA10", 35,  "P", "R4", "S1", "G1", 1.0,    30.0, 1),
		Dihedral_NA("DNA11", 36,  "P", "R4", "S1", "C1", 1.0,    30.0, 1),
		Dihedral_NA("DNA12", 37,  "P", "R4", "S1", "U1", 1.0,    30.0, 1),
		Dihedral_NA("DNA13", 38,  "C", "R4",  "P",  "O", 1.2,    5.0, 1),
		Dihedral_NA("DNA14", 39, "S1", "R4",  "P",  "O", 1.2,   160.0, 1),
		Dihedral_NA("DNA15", 40,  "O",  "C", "R4",  "P", 1.0,    90.0, 3),
		Dihedral_NA("DNA16", 41,  "O",  "C", "R4", "S1", 1.0,   150.0, 1),
		Dihedral_NA("DNA17", 42,  "P",  "O",  "C", "R4", 0.33,    0.0, 1),
		Dihedral_NA("DNA18", 43,  "P",  "O",  "C", "R4", 0.125, 180.0, 2),
		Dihedral_NA("DNA19", 44,  "P",  "O",  "C", "R4", 0.83,    0.0, 3),
		Dihedral_NA("DNA20", 45, "R4",  "P",  "O",  "C", 1.00,    0.0, 2)]

#old values
#       Dihedral_NA("DNA01", "R4", "S1", "G1", "G2", 1.0,  -120.0, 2),
#		Dihedral_NA("DNA02", "R4", "S1", "A1", "A2", 1.0,  -120.0, 2),
#		Dihedral_NA("DNA03", "R4", "S1", "A1", "A2", 0.2,   160.0, 1),    
#		Dihedral_NA("DNA04", "R4", "A1", "A2", "S1", 1.0,   165.0, 1),
#		Dihedral_NA("DNA05",  "C", "R4", "S1", "A1", 1.0,   150.0, 1),
#		Dihedral_NA("DNA06",  "C", "R4", "S1", "G1", 1.0,   150.0, 1),
#		Dihedral_NA("DNA07",  "C", "R4", "S1", "C1", 1.0,   150.0, 1),
#		Dihedral_NA("DNA08",  "C", "R4", "S1", "U1", 1.0,   150.0, 1),
#		Dihedral_NA("DNA09",  "P", "R4", "S1", "A1", 1.0,    30.0, 1),
#		Dihedral_NA("DNA10",  "P", "R4", "S1", "G1", 1.0,    30.0, 1),
#		Dihedral_NA("DNA11",  "P", "R4", "S1", "C1", 1.0,    30.0, 1),
#		Dihedral_NA("DNA12",  "P", "R4", "S1", "U1", 1.0,    30.0, 1),
#		Dihedral_NA("DNA13",  "C", "R4",  "P",  "O", 1.2,    5.0, 1),
#		Dihedral_NA("DNA14", "S1", "R4",  "P",  "O", 1.2,   160.0, 1),
#		Dihedral_NA("DNA15",  "O",  "C", "R4",  "P", 1.0,    90.0, 3),
#		Dihedral_NA("DNA16",  "O",  "C", "R4", "S1", 1.0,   150.0, 1),
#		Dihedral_NA("DNA17",  "P",  "O",  "C", "R4", 0.33,    0.0, 1),
#		Dihedral_NA("DNA18",  "P",  "O",  "C", "R4", 0.125, 180.0, 2),
#		Dihedral_NA("DNA19",  "P",  "O",  "C", "R4", 0.83,    0.0, 3),
#		Dihedral_NA("DNA20", "R4",  "P",  "O",  "C", 1.00,    0.0, 2)]