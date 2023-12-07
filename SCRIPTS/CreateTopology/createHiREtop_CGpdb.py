## @file createHiREtop_pdb.py
#
# Create topology from CG pdb file

import argparse
from os.path import exists
import PDBhandler
import FA2CGmapper
import topology
import ChemData

def create_start(coords):
    outf = open("start","w")
    for xyz in coords:
        line = "{:16.10f}".format(xyz[0]) + "{:16.10f}".format(xyz[1]) + "{:16.10f}".format(xyz[2]) + "\n"
        outf.write(line)
    outf.close()

def create_CGcoords(coordsbyres,nres):
    coords = list()
    for idres in range(nres):
        for xyz in coordsbyres[idres+1]:
            coords.append(xyz)
    return coords


#parse arguments
parser = argparse.ArgumentParser(prog = 'Create HiRE input from coarse grained pdb file',
                    description = 'This script creates a topology and start file from an coarse-grained pdb file',
                    epilog = 'Do not forget to get an up-to-date parameter file')

parser.add_argument('filename', help='Name of pdb file used to create HiRE input') 
args = parser.parse_args()
fname = args.filename
if not(exists(fname)):
    raise FileNotFoundError("File %s does not exist" % fname)
# get data from pdb
data = PDBhandler.parse_CGpdb(fname)
# unpack tuple
natom,nres,CG_labels,elements,CG_resrange,CG_resnames,coordsbyres,CG_termini = data

# get coords as array -> CG_coords
CG_coords = create_CGcoords(coordsbyres,nres)
# get CG particle masses
CG_masses = FA2CGmapper.assign_CGmasses(CG_labels)
# get CG particle charges and types
CG_charges, CG_partypes = FA2CGmapper.assign_par_properties(CG_labels)

# extension to study DNA and RNA: assign whether we have RNA or DNA!
moltype = ChemData.RNA_or_DNA(CG_resnames)

# translate all information into topology data dictionary
top_dict = topology.init_topology()
top_dict["nres"] = nres
top_dict["nparticles"] = natom
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
angles,angtype = ChemData.get_angles(top_dict["nmol"], CG_termini, CG_labels, moltype)
nangtypes,tk,teq,nangles,angles_top = ChemData.get_angleinfo(angles,angtype)
top_dict["nangles"] = nangles
top_dict["angs"] = angles_top
top_dict["nangtypes"] = nangtypes
top_dict["tk"] = tk
top_dict["teq"] = teq

#dihedral information
dihs,dihtypes = ChemData.get_dihs(top_dict["nmol"], CG_termini, CG_labels, moltype)
ntorstypes,pk,phi,ndihs,dihs_top,pn = ChemData.get_dihinfo(dihs,dihtypes)
top_dict["ndihs"] = ndihs
top_dict["dihs"] = dihs_top
top_dict["ndihstypes"] = ntorstypes
top_dict["pk"] = pk
top_dict["pn"] = pn
top_dict["phi"] = phi

topology.create_newtopology(top_dict)
create_start(CG_coords)