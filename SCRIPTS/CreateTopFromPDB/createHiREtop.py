import PDBhandler
import FA2CGmapper
import topology
import sys

data = PDBhandler.parse_pdb(sys.argv[1])
natom,nres,atomnames,elements,res,resnames,coordsbyres,termini = data
atommasses = FA2CGmapper.assign_atommass(atomnames)
resnames, CG_resnames = FA2CGmapper.assign_resnames(resnames)
CG_labels,CG_coords,CG_resrange = FA2CGmapper.get_CG_info(resnames,res,atomnames,atommasses,coordsbyres)
CG_masses = FA2CGmapper.assign_CGmasses(CG_labels)
CG_termini,CG_teri_flags = FA2CGmapper.get_CG_termini(res,termini,CG_resrange)
CG_charges, CG_partypes = FA2CGmapper.assign_par_properties(CG_labels)

top_dict = topology.init_topology()
top_dict["nres"] = nres
top_dict["nparticles"] = len(CG_labels)
top_dict["nmol"] = len(CG_termini)/2
top_dict["fl_mol"] = CG_termini
top_dict["particle_names"] = CG_labels
top_dict["res_names"] = CG_resnames
top_dict["res_start"] = [i[0] for i in CG_resrange]
top_dict["res_final"] = [i[1] for i in CG_resrange]
top_dict["particle_mass"] = CG_masses
top_dict["npartypes"] = len(list(set(CG_labels)))
top_dict["charges"] = CG_charges
top_dict["particle_type"] = CG_partypes


top_dict["nbonds"]
top_dict["nangles"]
top_dict["ndihs"]
top_dict["nbondtypes"]
top_dict["nangtypes"] 
top_dict["ndihstypes"] 

top_dict["rk"] = list()                 #bond spring constant
top_dict["req"] = list()                #equilibrium bond length
top_dict["tk"] = list()                 #angular spring constant
top_dict["teq"]  = list()               #equilibrium angle
top_dict["pk"] = list()                 #pk
top_dict["pn"] = list()                 #pn
top_dict["phi"]  = list()               #torsional phase
top_dict["bonds"] = list()              #bond information
top_dict["angs"]  = list()              #angle information
top_dict["dihs"] = list()               #torsional information
