class Atom_DNA:
    def __init__(self, idx, name, el, mass):
        self.idx = idx
        self.name = name
        self.el = el
        self.mass = mass

    def set_idx(self,newidx):
        self.idx = newidx

class Bond_DNA:
    def __init__(self, atom1, atom2, force, dist):
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

class Angle_DNA:
    def __init__(self, atom1, atom2, atom3, force, ang):
        self.atom1 = atom1
        self.atom2 = atom2
        self.atom3 = atom3
        self.force = force
        self.ang = ang
        
class Dihedral_DNA:
    def __init__(self, atom1, atom2, atom3, atom4, force, dih, nterm):
        self.atom1 = atom1
        self.atom2 = atom2
        self.atom3 = atom3
        self.atom4 = atom4
        self.force = force
        self.dih = dih
        self.nterm = nterm

DA = [	Atom_DNA(1,   "P",  "P",   36.970),
	Atom_DNA(2, "O5*",  "O",   16.000),
	Atom_DNA(3, "C5*",  "C",   12.010),
	Atom_DNA(4,  "CA", "R4",   20.000),
	Atom_DNA(5,  "CY", "R1",   12.000),
	Atom_DNA(6,  "A1", "A1",   67.000),
	Atom_DNA(7,  "A2", "A2",   67.000)]

DC = [	Atom_DNA(1,   "P",  "P",   36.970),
	Atom_DNA(2, "O5*",  "O",   16.000),
	Atom_DNA(3, "C5*",  "C",   12.010),
	Atom_DNA(4,  "CA", "R4",   20.000),
	Atom_DNA(5,  "CY", "R1",   12.000),
	Atom_DNA(6,  "C1", "C1",  130.000)]

DG = [	Atom_DNA(1,   "P",  "P",   36.970),
	Atom_DNA(2, "O5*",  "O",   16.000),
	Atom_DNA(3, "C5*",  "C",   12.010),
	Atom_DNA(4,  "CA", "R4",   20.000),
	Atom_DNA(5,  "CY", "R1",   12.000),
	Atom_DNA(6,  "G1", "G1",   75.000),
	Atom_DNA(7,  "G2", "G2",   75.000)]

DT = [	Atom_DNA(1,   "P",  "P",   36.970),
	Atom_DNA(2, "O5*",  "O",   16.000),
	Atom_DNA(3, "C5*",  "C",   12.010),
	Atom_DNA(4,  "CA", "R4",   20.000),
	Atom_DNA(5,  "CY", "R1",   12.000),
	Atom_DNA(6,  "T1", "T1",  131.000)]

DNA_atoms = [	Atom_DNA(0,   "P",  "P",   36.970),
		Atom_DNA(0, "O5*",  "O",   16.000),
		Atom_DNA(0, "C5*",  "C",   12.010),
		Atom_DNA(0,  "CA", "R4",   20.000),
		Atom_DNA(0,  "CY", "R1",   12.000),
		Atom_DNA(0,  "A1", "A1",   67.000),
		Atom_DNA(0,  "A2", "A2",   67.000),
		Atom_DNA(0,  "G1", "G1",   75.000),
		Atom_DNA(0,  "G2", "G2",   75.000),
		Atom_DNA(0,  "C1", "C1",  130.000),
		Atom_DNA(0,  "T1", "T1",  131.000),
		Atom_DNA(0,   "D",  "D", 1000.000)]

DNA_bonds = [	Bond_DNA("R4",  "P",  30, 3.950),
		Bond_DNA("R4", "R1", 200, 2.344),
		Bond_DNA("R1", "G1", 200, 2.622),
		Bond_DNA("R1", "A1", 200, 2.633),
		Bond_DNA("R1", "T1", 200, 3.062),
		Bond_DNA("R1", "C1", 200, 3.004),
		Bond_DNA("G1", "G2", 200, 2.450),
		Bond_DNA("A1", "A2", 200, 2.180),
		Bond_DNA( "C", "R4", 200, 1.520),
		Bond_DNA( "P",  "O", 200, 1.593),
		Bond_DNA( "O",  "C", 200, 1.430),
		Bond_DNA( "D",  "D",  40, 12.00),
		Bond_DNA("R4",  "D",  10, 12.00),
		Bond_DNA( "D",  "O",  10, 12.00),
		Bond_DNA( "D",  "P",  10, 12.00)]
		
DNA_angles = [	Angle_DNA("R4", "R1", "A1",  70, 130.0),
		Angle_DNA("R4", "R1", "T1",  70, 142.0),
		Angle_DNA("R4", "R1", "G1",  70, 131.0),
		Angle_DNA("R4", "R1", "C1",  70, 146.0),
		Angle_DNA("R1", "A1", "A2", 120, 116.0),
		Angle_DNA("R1", "G1", "G2", 120, 111.2),
		Angle_DNA( "P",  "O",  "C",  70, 120.0),
		Angle_DNA( "O",  "C", "R4",  70, 110.6),
		Angle_DNA( "C", "R4",  "P",  70, 110.0),
		Angle_DNA("R4",  "P",  "O",  50, 105.0),    
		Angle_DNA( "C", "R4", "R1",  70, 130.0),
		Angle_DNA("R1", "R4",  "P", 100,  85.0),
		Angle_DNA( "D",  "D",  "D",  80, 180.0)]
 
DNA_dihs = [	Dihedral_DNA("R4", "R1", "G1", "G2", 1.0,  -120.0, 2),
		Dihedral_DNA("R4", "R1", "A1", "A2", 1.0,  -120.0, 2),
		Dihedral_DNA("R4", "R1", "A1", "A2", 0.2,   160.0, 1),    
		Dihedral_DNA("R4", "A1", "A2", "R1", 1.0,   165.0, 1),
		Dihedral_DNA( "C", "R4", "R1", "A1", 1.0,   150.0, 1),
		Dihedral_DNA( "C", "R4", "R1", "G1", 1.0,   150.0, 1),
		Dihedral_DNA( "C", "R4", "R1", "C1", 1.0,   150.0, 1),
		Dihedral_DNA( "C", "R4", "R1", "U1", 1.0,   150.0, 1),
		Dihedral_DNA( "P", "R4", "R1", "A1", 1.0,    30.0, 1),
		Dihedral_DNA( "P", "R4", "R1", "G1", 1.0,    30.0, 1),
		Dihedral_DNA( "P", "R4", "R1", "C1", 1.0,    30.0, 1),
		Dihedral_DNA( "P", "R4", "R1", "U1", 1.0,    30.0, 1),
		Dihedral_DNA( "C", "R4",  "P",  "O", 1.2,    5.0, 1),
		Dihedral_DNA("R1", "R4",  "P",  "O", 1.2,   160.0, 1),
		Dihedral_DNA( "O",  "C", "R4",  "P", 1.0,    90.0, 3),
		Dihedral_DNA( "O",  "C", "R4", "R1", 1.0,   150.0, 1),
		Dihedral_DNA( "P",  "O",  "C", "R4", 0.33,    0.0, 1),
		Dihedral_DNA( "P",  "O",  "C", "R4", 0.125, 180.0, 2),
		Dihedral_DNA( "P",  "O",  "C", "R4", 0.83,    0.0, 3),
		Dihedral_DNA("R4",  "P",  "O",  "C", 1.00,    0.0, 2)] 

	
		 
