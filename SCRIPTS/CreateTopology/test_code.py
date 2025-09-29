import PDBhandler
import ChemData
import FA2CGmapper
import unittest

def natoms(pdbname):
    data = PDBhandler.parse_pdb(pdbname)
    return data[0]

def nres(pdbname):
    data = PDBhandler.parse_pdb(pdbname)
    return data[1]   

def resdic(pdbname):
    data = PDBhandler.parse_pdb(pdbname)
    return data[4]

def resnames(pdbname):
    data = PDBhandler.parse_pdb(pdbname)
    return data[5]

def termini(pdbname):
    data = PDBhandler.parse_pdb(pdbname)
    return data[7]

def simple_resnames(resnames):
    return FA2CGmapper.assign_resnames(resnames)[0]

def CG_resnames(resnames):
    return FA2CGmapper.assign_resnames(resnames)[1]

class Testclass(unittest.TestCase):
    #tests for PDBhandler
    def test_natoms(var):
        var.assertEqual(natoms("test_RNA.pdb"), 472)
    def test_nres(var):
        var.assertEqual(nres("test_RNA.pdb"), 15)
    def test_atom_termini(var):
        terids = [1, 157, 158, 472]
        var.assertListEqual(termini("test_RNA.pdb"), terids)
    def test_res_dic(var):
        res_dic = {1: [1, 28], 2: [29, 61], 3: [62, 92], 4: [93, 123], 5: [124, 157], 6: [158, 185], 7: [186, 219], 8: [220, 253], 9: [254, 283], 10: [284, 316], 11: [317, 346], 12: [347, 376], 13: [377, 407], 14: [408, 438], 15: [439, 472]}
        var.assertDictEqual(resdic("test_RNA.pdb"), res_dic)
    def test_res_names(var):
        res_names = ['U5', 'A', 'C', 'C', 'G3', 'C5', 'G', 'G', 'U', 'A', 'U', 'U', 'C', 'C', 'A3']
        var.assertListEqual(resnames("test_RNA.pdb"), res_names)
    #tests for FA2CGmapper
    def test_assign_mass(var):
        input_names = ["O3'", 'P', "C1'", "H1'", 'N1', 'H1', 'C6', 'O6']
        out_mass = [16.000,30.974,12.011,1.008,14.007,1.008,12.011,16.000] 
        var.assertListEqual(FA2CGmapper.assign_atommass(input_names),out_mass)
    def test_CG_masses(var):
        input_names = ["O5", "P", "O3", "R4"]
        var.assertListEqual(FA2CGmapper.assign_CGmasses(input_names),[16.000,36.970,16.000,20.000])
    def test_CG_termini(var):
        input_res = {1: [1,20], 2: [21,31], 3: [32,41]}
        input_termini = [1,31,32,41]
        input_CGres = {1:[1,6], 2:[7,13], 3:[14,21]}
        out_termini = [1,13,14,21]
        var.assertListEqual(FA2CGmapper.get_CG_termini(input_res, input_termini, input_CGres),out_termini)
    def test_simple_resnames(var):
        input_resnames = ['U5', 'A', 'C', 'C', 'G3']
        out_resnames = ['U', 'A', 'C', 'C', 'G']
        var.assertListEqual(simple_resnames(input_resnames),out_resnames)
    def test_simple_resnames(var):
        input_resnames = ['U5', 'A', 'C', 'C', 'G3']
        out_resnames = ['URAi', 'ADE', 'CYS', 'CYS', 'GUA']
        var.assertListEqual(CG_resnames(input_resnames),out_resnames)

if __name__ == '__main__':
    unittest.main()