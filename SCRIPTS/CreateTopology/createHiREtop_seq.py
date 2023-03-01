## @file createHiREtop_pdb.py
#
# Create topology from sequence
import argparse
from os.path import exists
import FA2CGmapper
import CGdatafromSeq
import topology
import ChemData

def read_seqlist(fname):
    seqlist = list()
    with open(fname, "r") as f:
        for line in f:
            seqlist += line.split()
    if len(seqlist) == 0:
        raise Exception("No sequence found in input file "+ fname)
    else:
        return seqlist

#parse arguments
parser = argparse.ArgumentParser(prog = 'Create HiRE input from sequence data',
                    description = 'This script creates a topology from sequence',
                    epilog = 'You can get a start file from a CG pdb file')

parser.add_argument('filename', help='Name of sequence file used to create HiRE input') 
args = parser.parse_args()
fname = args.filename
if not(exists(fname)):
    raise FileNotFoundError("File %s does not exist" % fname)
# read sequence
seqlist = read_seqlist(fname)
# get list of atoms
CG_labels = CGdatafromSeq.get_atoms_from_seq(seqlist)
# get res names used for topology
CG_resnames = seqlist
# use call below to get old style topology
# CG_resnames = CGdatafromSeq.translate_CGresnames(seqlist)
# get residue ranges
CG_resrange = CGdatafromSeq.get_resrange(seqlist)
# get masses
CG_masses = FA2CGmapper.assign_CGmasses(CG_labels)
# get CG molecule termini
CGres_termini = CGdatafromSeq.find_termini_seqlist(seqlist)
CG_termini = CGdatafromSeq.get_CGterminis(CG_resrange, CGres_termini)
# finally, get CG particle charges and types
CG_charges, CG_partypes = FA2CGmapper.assign_par_properties(CG_labels)

nres = len(seqlist)

# translate all information into topology data dictionary
top_dict = topology.init_topology()
top_dict["nres"] = nres
top_dict["nparticles"] = len(CG_labels)
top_dict["nmol"] = int(len(CG_termini)/2)
top_dict["fl_mol"] = CG_termini
top_dict["particle_names"] = CG_labels
top_dict["res_names"] = CG_resnames
top_dict["res_start"] = [CG_resrange[i][0] for i in CG_resrange.keys()]
top_dict["res_final"] = [CG_resrange[i][1] for i in CG_resrange.keys()]
top_dict["particle_mass"] = CG_masses
top_dict["npartypes"] = len(list(set(CG_labels)))
top_dict["charges"] = CG_charges
top_dict["particle_type"] = CG_partypes

#bonding information
bonds,bondtype = ChemData.get_bonds(top_dict["nmol"], CG_termini, CG_labels)
nbondtypes,rk,req,nbonds,bonds_top = ChemData.get_bondinfo(bonds,bondtype)
top_dict["nbonds"] = nbonds
top_dict["bonds"] = bonds_top
top_dict["nbondtypes"] = nbondtypes
top_dict["rk"] = rk
top_dict["req"] = req

#angle information
angles,angtype = ChemData.get_angles(top_dict["nmol"], CG_termini, CG_labels)
nangtypes,tk,teq,nangles,angles_top = ChemData.get_angleinfo(angles,angtype)
top_dict["nangles"] = nangles
top_dict["angs"] = angles_top
top_dict["nangtypes"] = nangtypes
top_dict["tk"] = tk
top_dict["teq"] = teq

#dihedral information
dihs,dihtypes = ChemData.get_dihs(top_dict["nmol"], CG_termini, CG_labels)
ntorstypes,pk,phi,ndihs,dihs_top,pn = ChemData.get_dihinfo(dihs,dihtypes)
top_dict["ndihs"] = ndihs
top_dict["dihs"] = dihs_top
top_dict["ndihstypes"] = ntorstypes
top_dict["pk"] = pk
top_dict["pn"] = pn
top_dict["phi"] = phi

topology.create_newtopology(top_dict)