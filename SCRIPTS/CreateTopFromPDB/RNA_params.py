class Atom_RNA:
    def __init__(self, idx, name, el, mass):
    	self.idx = idx
    	self.name = name
    	self.el = el
    	self.mass = mass

    def set_idx(self,newidx):
    	self.idx = newidx

class Bond_RNA:
    def __init__(self, atom1, atom2, force, dist):
        self.atom1 = atom1
        self.atom2 = atom2
        self.force = force
        self.dist = dist
    
    def iscorrectbond(self, at1, at2):
    	if ((at1 is self.atom1) and (at2 is self.atom2)):
    	    return True
    	elif ((at12 is self.atom1) and (at1 is self.atom2)):
    	    return True
    	else:
    	    return False

class Angle_RNA:
    def __init__(self, atom1, atom2, atom3, force, ang):
        self.atom1 = atom1
        self.atom2 = atom2
        self.atom3 = atom3
        self.force = force
        self.ang = ang
        
class Dihedral_RNA:
    def __init__(self, atom1, atom2, atom3, atom4, force, dih, nterm):
        self.atom1 = atom1
        self.atom2 = atom2
        self.atom3 = atom3
        self.atom4 = atom4
        self.force = force
        self.dih = dih
        self.nterm = nterm

RA = [	Atom_RNA(1,   "P",  "P",   36.970),
	Atom_RNA(2, "O5*",  "O",   16.000),
	Atom_RNA(3, "C5*",  "C",   12.010),
	Atom_RNA(4,  "CA", "R4",   20.000),
	Atom_RNA(5,  "CY", "R1",   12.000),
	Atom_RNA(6,  "A1", "A1",   67.000),
	Atom_RNA(7,  "A2", "A2",   67.000)]

RC = [	Atom_RNA(1,   "P",  "P",   36.970),
	Atom_RNA(2, "O5*",  "O",   16.000),
	Atom_RNA(3, "C5*",  "C",   12.010),
	Atom_RNA(4,  "CA", "R4",   20.000),
	Atom_RNA(5,  "CY", "R1",   12.000),
	Atom_RNA(6,  "C1", "C1",  130.000)]

RG = [	Atom_RNA(1,   "P",  "P",   36.970),
	Atom_RNA(2, "O5*",  "O",   16.000),
	Atom_RNA(3, "C5*",  "C",   12.010),
	Atom_RNA(4,  "CA", "R4",   20.000),
	Atom_RNA(5,  "CY", "R1",   12.000),
	Atom_RNA(6,  "G1", "G1",   75.000),
	Atom_RNA(7,  "G2", "G2",   75.000)]

RU = [	Atom_RNA(1,   "P",  "P",   36.970),
	Atom_RNA(2, "O5*",  "O",   16.000),
	Atom_RNA(3, "C5*",  "C",   12.010),
	Atom_RNA(4,  "CA", "R4",   20.000),
	Atom_RNA(5,  "CY", "R1",   12.000),
	Atom_RNA(6,  "U1", "U1",  131.000)]

RNA_atoms = [	Atom_RNA(0,   "P",  "P",   36.970),
		Atom_RNA(0, "O5*",  "O",   16.000),
		Atom_RNA(0, "C5*",  "C",   12.010),
		Atom_RNA(0,  "CA", "R4",   20.000),
		Atom_RNA(0,  "CY", "R1",   12.000),
		Atom_RNA(0,  "A1", "A1",   67.000),
		Atom_RNA(0,  "A2", "A2",   67.000),
		Atom_RNA(0,  "G1", "G1",   75.000),
		Atom_RNA(0,  "G2", "G2",   75.000),
		Atom_RNA(0,  "C1", "C1",  130.000),
		Atom_RNA(0,  "U1", "U1",  131.000),
		Atom_RNA(0,   "D",  "D", 1000.000)]

RNA_bonds = [	Bond_RNA("R4",  "P",  30, 3.800),
		Bond_RNA("R4", "R1", 200, 2.344),
		Bond_RNA("R1", "G1", 200, 2.622),
		Bond_RNA("R1", "A1", 200, 2.633),
		Bond_RNA("R1", "U1", 200, 3.062),
		Bond_RNA("R1", "C1", 200, 3.004),
		Bond_RNA("G1", "G2", 200, 2.450),
		Bond_RNA("A1", "A2", 200, 2.180),
		Bond_RNA( "C", "R4", 200, 1.520),
		Bond_RNA( "P",  "O", 200, 1.593),
		Bond_RNA( "O",  "C", 200, 1.430),
		Bond_RNA( "D",  "D",  40, 12.00),
		Bond_RNA("R4",  "D",  10, 12.00),
		Bond_RNA( "D",  "O",  10, 12.00),
		Bond_RNA( "D",  "P",  10, 12.00)]
		
RNA_angles = [	Angle_RNA("R4", "R1", "A1",  70, 123.6),
		Angle_RNA("R4", "R1", "U1",  70, 132.5),
		Angle_RNA("R4", "R1", "G1",  70, 123.5),
		Angle_RNA("R4", "R1", "C1",  70, 131.1),
		Angle_RNA("R1", "A1", "A2", 120, 116.7),
		Angle_RNA("R1", "G1", "G2", 120, 111.2),
		Angle_RNA( "P",  "O",  "C",  70, 122.9),
		Angle_RNA( "O",  "C", "R4",  70, 110.6),
		Angle_RNA( "C", "R4",  "P",  70,  98.0),
		Angle_RNA("R4",  "P",  "O",  50, 110.0),    
		Angle_RNA( "C", "R4", "R1",  70, 135.5),
		Angle_RNA("R1", "R4",  "P", 100,  98.0),
		Angle_RNA( "D",  "D",  "D",  80, 180.0)]
 
RNA_dihs = [	Dihedral_RNA("R4", "R1", "G1", "G2", 1.0,     -20, 1),
		Dihedral_RNA("R4", "R1", "G1", "G2", 0.2,    -160, 1),
		Dihedral_RNA("R4", "R1", "A1", "A2", 1.0,     -20, 1),
		Dihedral_RNA("R4", "R1", "A1", "A2", 0.2,  -160.0, 1),    
		Dihedral_RNA("R4", "A1", "A2", "R1", 1.0,  -165.0, 1),
		Dihedral_RNA("R4", "G1", "G2", "R1", 1.0,  -165.0, 1),
		Dihedral_RNA( "C", "R4", "R1", "A1", 1.0,  -150.0, 1),
		Dihedral_RNA( "C", "R4", "R1", "A1", 0.2,   150.0, 1),
		Dihedral_RNA( "C", "R4", "R1", "G1", 1.0,  -150.0, 1),
		Dihedral_RNA( "C", "R4", "R1", "C1", 1.0,  -150.0, 1),
		Dihedral_RNA( "C", "R4", "R1", "C1", 0.2,   150.0, 1),
		Dihedral_RNA( "C", "R4", "R1", "U1", 1.0,  -150.0, 1),
		Dihedral_RNA( "C", "R4", "R1", "U1", 0.2,   150.0, 1),
		Dihedral_RNA( "P", "R4", "R1", "A1", 1.0,   100.0, 1),
		Dihedral_RNA( "P", "R4", "R1", "G1", 1.0,   100.0, 1),
		Dihedral_RNA( "P", "R4", "R1", "C1", 1.0,   100.0, 1),
		Dihedral_RNA( "P", "R4", "R1", "U1", 1.0,   100.0, 1),
		Dihedral_RNA( "C", "R4",  "P",  "O", 1.2,    20.0, 1),
		Dihedral_RNA("R1", "R4",  "P",  "O", 1.2,   160.0, 1),
		Dihedral_RNA( "O",  "C", "R4",  "P", 1.0,    90.0, 3),
		Dihedral_RNA( "O",  "C", "R4", "R1", 1.0,   150.0, 1),
		Dihedral_RNA( "P",  "O",  "C", "R4", 0.33,    0.0, 1),
		Dihedral_RNA( "P",  "O",  "C", "R4", 0.125, 180.0, 2),
		Dihedral_RNA( "P",  "O",  "C", "R4", 0.83,    0.0, 3),
		Dihedral_RNA("R4",  "P",  "O",  "C", 1.00,    0.0, 2)] 

	
		 
