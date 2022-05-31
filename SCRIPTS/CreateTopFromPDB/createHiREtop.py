import PDBhandler
import FA2CGmapper
import topology
import sys

data = PDBhandler.parse_pdb(sys.argv[1])
natom,nres,atomnames,elements,res,resnames,coordsbyres,termini = data
atommasses = FA2CGmapper.assign_atommass(atomnames)
CG_labels,CG_coords,CG_resrange = FA2CGmapper.get_CG_info(resnames,res,atomnames,atommasses,coordsbyres)
CG_masses = FA2CGmapper.assign_CGmasses(CG_labels)
CG_termini = FA2CGmapper.get_CG_termini(res,termini,CG_resrange)

top_dict = topology.init_topology()
top_dict["nres"] = nres
top_dict["nparticles"] = len(CG_labels)
top_dict["nmol"] = len(CG_termini)/2
top_dict["fl_mol"] = CG_termini
top_dict["particle_names"] = CG_labels
top_dict["res_names"] = list()