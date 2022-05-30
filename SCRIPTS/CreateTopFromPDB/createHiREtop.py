import PDBhandler
import FA2CGmapper
import sys

data = PDBhandler.parse_pdb(sys.argv[1])
natom,nres,atomnames,elements,res,resnames,coordsbyres,termini = data
CG_labels,CG_coords,CG_resrange = FA2CGmapper.get_CG_info(resnames,res,atomnames,atommasses,coordsbyres)

