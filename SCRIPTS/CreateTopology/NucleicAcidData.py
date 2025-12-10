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
    
class QAngle_NA:
    def __init__(self, name, typeid, atom1, atom2, atom3, k, ref, a, b, c, d, e):
        self.name = name
        self.typeid = typeid
        self.atom1 = atom1
        self.atom2 = atom2
        self.atom3 = atom3
        self.k = k
        self.ref = ref        
        self.a = a
        self.b = b
        self.c = c
        self.d = d
        self.e = e

    def iscorrectangle(self, at1, at2, at3):
        if ((at1 is self.atom1) and (at2 is self.atom2) and (at3 is self.atom3)):
            return True
        else:
            return False
        
    def get_typeid(self):
        return self.typeid        

        
class Dihedral_NA:
    def __init__(self, name, typeid, atom1, atom2, atom3, atom4, force, dih, nterm, yoff):
        self.name = name
        self.typeid = typeid
        self.atom1 = atom1
        self.atom2 = atom2
        self.atom3 = atom3
        self.atom4 = atom4
        self.force = force
        self.dih = dih
        self.nterm = nterm
        self.yoff = yoff
    
    def iscorrectdih(self, at1, at2, at3, at4):
        if ((at1 is self.atom1) and (at2 is self.atom2) and (at3 is self.atom3) and (at4 is self.atom4)):
            return True
        else:
            return False  

    def get_typeid(self):
        return self.typeid	      

RA = [	Atom_NA(1,   "P",  "P",   36.970),
	Atom_NA(2, "O5*", "O5",   16.000),
	Atom_NA(3,  "CA", "R4",   20.000),
    Atom_NA(4, "O3*", "O3",   16.000),
	Atom_NA(5,  "CY", "R1",   12.000),
	Atom_NA(6,  "A1", "A1",   67.000),
	Atom_NA(7,  "A2", "A2",   67.000)]

DA = [	Atom_NA(1,   "P",  "P",   36.970),
	Atom_NA(2, "O5*",  "O",   16.000),
	Atom_NA(3,  "CA", "S4",   20.000),
    Atom_NA(4, "O3*", "O3",   16.000),
	Atom_NA(5,  "CY", "S1",   12.000),
	Atom_NA(6,  "A1", "A1",   67.000),
	Atom_NA(7,  "A2", "A2",   67.000)]

RC = [	Atom_NA(1,   "P",  "P",   36.970),
	Atom_NA(2, "O5*",  "O",   16.000),
	Atom_NA(3,  "CA", "R4",   20.000),
    Atom_NA(4, "O3*", "O3",   16.000),
	Atom_NA(5,  "CY", "R1",   12.000),
	Atom_NA(6,  "C1", "C1",  130.000)]

DC = [	Atom_NA(1,   "P",  "P",   36.970),
	Atom_NA(2, "O5*",  "O",   16.000),
	Atom_NA(3,  "CA", "S4",   20.000),
    Atom_NA(4, "O3*", "O3",   16.000),
	Atom_NA(5,  "CY", "S1",   12.000),
	Atom_NA(6,  "C1", "C1",  130.000)]

RG = [	Atom_NA(1,   "P",  "P",   36.970),
	Atom_NA(2, "O5*",  "O",   16.000),
	Atom_NA(3,  "CA", "R4",   20.000),
    Atom_NA(4, "O3*", "O3",   16.000),
	Atom_NA(5,  "CY", "R1",   12.000),
	Atom_NA(6,  "G1", "G1",   75.000),
	Atom_NA(7,  "G2", "G2",   75.000)]

DG = [	Atom_NA(1,   "P",  "P",   36.970),
	Atom_NA(2, "O5*",  "O",   16.000),
	Atom_NA(3,  "CA", "S4",   20.000),
    Atom_NA(4, "O3*", "O3",   16.000),
	Atom_NA(5,  "CY", "S1",   12.000),
	Atom_NA(6,  "G1", "G1",   75.000),
	Atom_NA(7,  "G2", "G2",   75.000)]

DT = [	Atom_NA(1,   "P",  "P",   36.970),
	Atom_NA(2, "O5*",  "O",   16.000),
	Atom_NA(3,  "CA", "S4",   20.000),
    Atom_NA(4, "O3*", "O3",   16.000),
	Atom_NA(5,  "CY", "S1",   12.000),
	Atom_NA(6,  "T1", "T1",  131.000)]

RU = [	Atom_NA(1,   "P",  "P",   36.970),
	Atom_NA(2, "O5*",  "O",   16.000),
	Atom_NA(3,  "CA", "R4",   20.000),
    Atom_NA(4, "O3*", "O3",   16.000),
	Atom_NA(5,  "CY", "R1",   12.000),
	Atom_NA(6,  "U1", "U1",  131.000)]

NA_atoms = [	Atom_NA(0,   "P",  "P",   36.970),
		Atom_NA(0, "O5*", "O5",   16.000),
		Atom_NA(0,  "CA", "R4",   20.000),  # RNA sugar
        Atom_NA(0,  "CA", "S4",   20.000),  # DNA sugar
		Atom_NA(0,  "CY", "R1",   12.000),  # RNA sugar
 		Atom_NA(0,  "CY", "S1",   12.000),  # DNA sugar
		Atom_NA(0, "O3*", "O3",   16.000),              
		Atom_NA(0,  "A1", "A1",   67.000),
		Atom_NA(0,  "A2", "A2",   67.000),
		Atom_NA(0,  "G1", "G1",   75.000),
		Atom_NA(0,  "G2", "G2",   75.000),
		Atom_NA(0,  "C1", "C1",  130.000),
        Atom_NA(0,  "T1", "T1",  131.000),
		Atom_NA(0,  "U1", "U1",  131.000),
		Atom_NA(0,   "D",  "D", 1000.000)]

NA_bonds = [	
        Bond_NA( 1, "R4", "O3", 200, 2.4647), #RNA
 		Bond_NA( 2, "S4", "O3", 200, 2.4647), #DNA
		Bond_NA( 3, "R4", "R1", 200, 2.3242), #RNA
		Bond_NA( 4, "S4", "S1", 200, 2.3242), #DNA   
		Bond_NA( 5, "R1", "G1", 200, 2.6601), #RNA
		Bond_NA( 6, "R1", "A1", 200, 2.6538), #RNA
		Bond_NA( 7, "R1", "U1", 200, 3.0943), #RNA
		Bond_NA( 8, "R1", "C1", 200, 3.0310), #RNA
		Bond_NA( 9, "S1", "G1", 200, 2.622), #DNA
		Bond_NA(10, "S1", "A1", 200, 2.633), #DNA
		Bond_NA(11, "S1", "T1", 200, 3.062), #DNA
		Bond_NA(12, "S1", "C1", 200, 3.004), #DNA   
		Bond_NA(13, "G1", "G2", 200, 2.4859),
		Bond_NA(14, "A1", "A2", 200, 2.2290),
		Bond_NA(15, "O5", "R4", 200, 2.4591), #RNA
        Bond_NA(16, "O5", "S4", 200, 2.4591), #DNA
		Bond_NA(17, "P",  "O3", 200, 1.6089),
		Bond_NA(18, "O5",  "P", 200, 1.6153),
		Bond_NA(19,  "D",  "D",  40, 12.00),
		Bond_NA(20, "R4",  "D",  10, 12.00),
        Bond_NA(21, "S4",  "D",  10, 12.00),
		Bond_NA(22,  "D", "O3",  10, 12.00),
        Bond_NA(23,  "D", "O5",  10, 12.00),
		Bond_NA(24,  "D",  "P",  10, 12.00)]


NA_angles = [	
        Angle_NA("RNA01",  1, "R4", "R1", "A1",  70, 123.6),
		Angle_NA("RNA02",  2, "R4", "R1", "U1",  70, 132.5),
		Angle_NA("RNA03",  3, "R4", "R1", "G1",  70, 120.2), #new value
		Angle_NA("RNA04",  4, "R4", "R1", "C1",  70, 131.1), 
		Angle_NA("RNA05",  5, "R1", "A1", "A2", 120, 115.1), #new value
		Angle_NA("RNA06",  6, "R1", "G1", "G2", 120, 110.8), #new value
		Angle_NA("RNA07",  7,  "P", "O5", "R4",  70, 156.9), #new value
		Angle_NA("RNA08",  8, "O5", "R4", "O3",  70, 110.6),
		Angle_NA("RNA09",  9, "R4", "O3",  "P",  70,  98.0),
		Angle_NA("RNA10", 10, "O3",  "P", "O5",  50, 104.5), #new value
		Angle_NA("RNA11", 11, "O5", "R4", "R1",  70, 135.5),
		Angle_NA("RNA12", 12, "R1", "R4", "O3", 100,  98.0),
		Angle_NA("RNA13", 13,  "D",  "D",  "D",  80, 180.0)]

NA_qangles = [
    
]

NA_dihs = [	
    	Dihedral_NA("RNA01",  1, "R4", "R1", "A1", "A2", 30.0, 6.0, 1, 26.2),
        Dihedral_NA("RNA02",  2, "R4", "R1", "A1", "A2",  3.5, 2.0, 5, 26.2),
    	Dihedral_NA("RNA03",  3, "R4", "R1", "G1", "G2", 35.7, 0.0, 1, 23.6),
        Dihedral_NA("RNA04",  4, "R4", "R1", "G1", "G2",  6.0, 6.0, 4, 23.6),        
		Dihedral_NA("RNA05",  5, "R4", "A1", "A2", "R1", 60.0, 3.1, 1, 17.0),
		Dihedral_NA("RNA06",  6, "R4", "G1", "G2", "R1", 60.0, 3.1, 1, 20.0),
		Dihedral_NA("RNA07",  7, "O5", "R4", "R1", "A1", 33.0, 3.1, 1, 24.7),
		Dihedral_NA("RNA08",  8, "O5", "R4", "R1", "A1",  4.4, 4.6, 6, 24.7),        
		Dihedral_NA("RNA09",  9, "O5", "R4", "R1", "A1",  5.0, 4.6, 7, 24.7),
		Dihedral_NA("RNA10", 10, "O5", "R4", "R1", "C1", 32.0, 3.1, 1, 27.0),
		Dihedral_NA("RNA11", 11, "O5", "R4", "R1", "C1",  8.7, 4.1, 6, 27.0),        
		Dihedral_NA("RNA12", 12, "O5", "R4", "R1", "C1",  2.7, 5.8, 7, 27.0),
        Dihedral_NA("RNA13", 13, "O5", "R4", "R1", "G1", 33.0, 3.1, 1, 24.7),
		Dihedral_NA("RNA14", 14, "O5", "R4", "R1", "G1",  4.4, 4.6, 6, 24.7),        
		Dihedral_NA("RNA15", 15, "O5", "R4", "R1", "G1",  5.0, 4.6, 7, 24.7),
        Dihedral_NA("RNA16", 16, "O5", "R4", "R1", "U1", 33.8, 3.1, 1, 29.0),
		Dihedral_NA("RNA17", 17, "O5", "R4", "R1", "U1",  6.4, 4.0, 6, 29.0),        
		Dihedral_NA("RNA18", 18, "O5", "R4", "R1", "U1",  3.0, 4.7, 7, 29.0),   
        Dihedral_NA("RNA19", 19, "O3", "R4", "R1", "A1", 14.4, 0.7, 1, 20.0),
		Dihedral_NA("RNA20", 20, "O3", "R4", "R1", "A1", 14.7, 0.6, 1, 20.0),        
		Dihedral_NA("RNA21", 21, "O3", "R4", "R1", "A1",  5.9, 3.2, 4, 20.0),
		Dihedral_NA("RNA22", 22, "O3", "R4", "R1", "C1", 21.5, 1.1, 1, 22.0),
		Dihedral_NA("RNA23", 23, "O3", "R4", "R1", "C1", 16.0, 0.6, 1, 22.0),        
		Dihedral_NA("RNA24", 24, "O3", "R4", "R1", "C1",  9.8, 3.2, 4, 22.0),
        Dihedral_NA("RNA25", 25, "O3", "R4", "R1", "G1", 15.2, 1.4, 1, 21.0),
		Dihedral_NA("RNA26", 26, "O3", "R4", "R1", "G1", 17.0, 0.6, 1, 21.0),        
		Dihedral_NA("RNA27", 27, "O3", "R4", "R1", "G1",  7.9, 3.2, 4, 21.0),
        Dihedral_NA("RNA28", 28, "O3", "R4", "R1", "U1", 15.2, 1.1, 1, 21.2),
		Dihedral_NA("RNA29", 29, "O3", "R4", "R1", "U1", 16.7, 0.6, 1, 21.2),        
		Dihedral_NA("RNA30", 30, "O3", "R4", "R1", "U1",  9.3, 3.2, 4, 21.2),      
        Dihedral_NA("RNA31", 31, "O5", "R4", "O3",  "P",  8.0, 3.6, 1,  0.0), 
        Dihedral_NA("RNA32", 32, "O5", "R4", "O3",  "P",  4.0, 1.5, 3,  0.0),        
        Dihedral_NA("RNA33", 33, "O5", "R4", "O3",  "P",  2.0, 2.0, 4,  0.0),
        Dihedral_NA("RNA34", 34, "R1", "R4", "O3",  "P",  5.8, 5.5, 1,  5.7),        
        Dihedral_NA("RNA35", 35, "R1", "R4", "O3",  "P",  3.0, 1.2, 3,  5.7), 
        Dihedral_NA("RNA36", 36, "R1", "R4", "O3",  "P",  2.0, 3.3, 4,  5.7), 
        Dihedral_NA("RNA37", 37, "P" , "O5", "R4", "O3",  4.7, 0.0, 1,  0.0),
        Dihedral_NA("RNA38", 38, "P" , "O5", "R4", "O3",  5.0, 5.6, 2,  0.0),
        Dihedral_NA("RNA39", 39, "P" , "O5", "R4", "R1",  2.0, 6.1, 1, 12.8),
        Dihedral_NA("RNA40", 40, "P" , "O5", "R4", "R1",  2.8, 1.4, 2, 12.8),        
        Dihedral_NA("RNA41", 41, "O3", "P" , "O5", "R4",  1.4, 3.4, 1,  0.0),
        Dihedral_NA("RNA42", 42, "O3", "P" , "O5", "R4",  3.5, 6.9, 2,  0.0),
        Dihedral_NA("RNA43", 43, "O3", "P" , "O5", "R4",  1.6, 1.0, 3,  0.0),
        Dihedral_NA("RNA44", 44, "R4", "O3", "P" , "O5",  3.3, 2.2, 1,  0.0),    
        Dihedral_NA("RNA45", 45, "R4", "O3", "P" , "O5",  3.5, 5.3, 2,  0.0)]

