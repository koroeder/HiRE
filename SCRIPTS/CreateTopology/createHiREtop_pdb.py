## @file createHiREtop_pdb.py
#
# Create topology from FA pdb file

import PDBhandler
import FA2CGmapper
import topology
import ChemData
import sys

def create_start(coords):
    outf = open("start","w")
    for xyz in coords:
        line = "{:16.10f}".format(xyz[0]) + "{:16.10f}".format(xyz[1]) + "{:16.10f}".format(xyz[2]) + "\n"
        outf.write(line)
    outf.close()

# get data from pdb
data = PDBhandler.parse_pdb(sys.argv[1])
# unpack tuple
natom,nres,atomnames,elements,res,resnames,coordsbyres,termini = data
# get masses for all atoms
atommasses = FA2CGmapper.assign_atommass(atomnames)
# get CG res names and simplified res names
resnames, CG_resnames = FA2CGmapper.assign_resnames(resnames)
# obtain CG particle labels and coordinates, and map CG residue start and final particles
CG_labels,CG_coords,CG_resrange = FA2CGmapper.get_CG_info(resnames,res,atomnames,atommasses,coordsbyres)
# get CG particle masses
CG_masses = FA2CGmapper.assign_CGmasses(CG_labels)
# get CG molecule termini
CG_termini = FA2CGmapper.get_CG_termini(res,termini,CG_resrange)
# finally, get CG particle charges and types
CG_charges, CG_partypes = FA2CGmapper.assign_par_properties(CG_labels)

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
create_start(CG_coords)