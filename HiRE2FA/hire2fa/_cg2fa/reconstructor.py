import json
import numpy as np
from pathlib import Path

import hire2fa as h2f

# //////////////////////////////////////////////////////////////////////////////
class Reconstructor:
    def __init__(self, path_pdb_cg: Path, path_model: Path):
        pdb = h2f.PDBTable.read_pdb(path_pdb_cg)
        self.chains = h2f.Residue.get_pdb_chains(pdb, init_with_cg_data = True)
        self.residues = [r for chain in self.chains for r in chain]

        model_raw: dict[str, list] = json.loads(path_model.read_text())

        self.model: dict[ str, dict[str, list[list[float]]] ] = {
            "_" : {}, # generic atoms (backbone and sugar)
            "U" : {}, "C" : {}, "A" : {}, "G" : {},
        }

        for name,lst in model_raw.items():
            translated = h2f.Mapping.MAP_STATSNAME_TO_STANDARDNAME[name]
            if   name.startswith("bu"): key = 'U'
            elif name.startswith("bc"): key = 'C'
            elif name.startswith("ba"): key = 'A'
            elif name.startswith("bg"): key = 'G'
            else: key = '_'
            self.model[key][translated] = lst


    # ------------------------------------------------------------------------------
    def reconstruct(self) -> None:
        self._apply_direct_mappings() # P, O3', O5', C1', C4'
        self._reconstruct_phosphate() # OP1, OP2
        self._reconstruct_backbone_and_sugar() # C3', C5', C2', O4', O2'
        self._reconstruct_nitro_base() #  N1, C2, C6 (pyr) / N9, C8, C4 (pur)
        self._drop_temp_particles() # remove the temporary centroids used for geometry calculations


    # --------------------------------------------------------------------------
    def export_reconstructed(self, path_pdb_fa: Path) -> None:
        out = h2f.PDBTable.join(r.fa_data for r in self.residues)
        out.normalize_indices()
        out.write_pdb(path_pdb_fa)


    # --------------------------------------------------------------------------
    def _apply_direct_mappings(self) -> None:
        """Starts populating `fa_data` with some direct mappings from `cg_data`."""
        centroids = [h2f.Mapping.name_centroid(1), h2f.Mapping.name_centroid(2)]
        for residue in self.residues:
            for i,(cg_name, fa_name) in enumerate(h2f.Mapping.iter_direct_mapping()):
                if cg_name == "P" and residue.is_5terminus: continue
                if fa_name in centroids and cg_name[0] != residue.resname[0]: continue

                particle = residue.cg_data.copy_filtered(atomname = cg_name)

                particle.sections["atomid"] = (str(i),)
                particle.sections["atomname"] = (fa_name,) # particle sections are guaranteed to have a single element (thanks to _assert_cg_residue)
                particle.sections["segid"] = particle.sections["chainid"]
                particle.sections["element"] = (fa_name[0],)
                residue.fa_data.append_table(particle)


    # ------------------------------------------------------------------------------
    def _reconstruct_phosphate(self) -> None:
        """
        reconstruct the phosphate oxygens. They have only 1 variant, so the reconstruction is well behaved.
        the O3 atoms are shifted in a copy of the residues to facilitate geometry calculations
        """
        chains_shift = [
            [residue.copy() for residue in chain]
            for chain in self.chains
        ]
        h2f.Residue.apply_o3_shift(chains_shift, do_cg = False)
        residues_shift = [r for chain in chains_shift for r in chain]

        names_phosph = ("O3'", "P", "O5'")
        for r_orig, r_shift in zip(self.residues, residues_shift):
            geo_ref = r_shift.init_fa_geometries(*names_phosph)
            pos_op1 = self._get_new_pos(geo_ref, "OP1", 0)
            pos_op2 = self._get_new_pos(geo_ref, "OP2", 0)
            self._reconstruct_fa_particle(r_orig, "OP1", pos_op1)
            self._reconstruct_fa_particle(r_orig, "OP2", pos_op2)


    # ------------------------------------------------------------------------------
    def _reconstruct_backbone_and_sugar(self) -> None:
        def get_coords(r: h2f.Residue, atomname: str) -> np.ndarray | None:
            return r.fa_data.copy_filtered(atomname = atomname).get_coords_array()[0]

        def infer_posits(refs: list[np.ndarray], variants: list[int]):
            pos_o5, pos_c1, pos_c4 = refs
            v_c5, v_o4, v_c3, v_c2, v_o2 = variants

            geo_c5 = h2f.Geometry(pos_c1, pos_c4, pos_o5)
            pos_c5 = self._get_new_pos(geo_c5, "C5'", v_c5)

            geo_o4 = h2f.Geometry(pos_c5, pos_c4, pos_c1)
            pos_o4 = self._get_new_pos(geo_o4, "O4'", v_o4)

            geo_c3 = h2f.Geometry(pos_o5, pos_c4, pos_c1)
            pos_c3 = self._get_new_pos(geo_c3, "C3'", v_c3)

            geo_o2 = h2f.Geometry(pos_c4, pos_o4, pos_c1)
            pos_c2 = self._get_new_pos(geo_o2, "C2'", v_c2)

            geo_o2 = h2f.Geometry(pos_c4, pos_c1, pos_c2)
            pos_o2 = self._get_new_pos(geo_o2, "O2'", v_o2)

            return pos_c5, pos_o4, pos_c3, pos_c2, pos_o2

        n_targets = 5
        nvars_per_target = 2
        combinations = nvars_per_target ** n_targets
        xyz_coords = 3

        all_variants = [
            [
                (i // (nvars_per_target ** j)) % nvars_per_target
                for j in range(n_targets)
            ]
            for i in range(combinations)
        ]

        avg_c1_c2_dist = 1.54
        avg_c2_c3_dist = 1.53
        avg_c3_c4_dist = 1.52
        avg_c4_o4_dist = 1.45
        avg_o4_c1_dist = 1.42
        avg_c3_o3_dist = 1.4
        avg_o2_o3_dist = 2.7
        avg_o2_c3_dist = 2.44

        for residue in self.residues:
            pos_o3_input = get_coords(residue, "O3'")
            pos_o5_input = get_coords(residue, "O5'")
            pos_c1_input = get_coords(residue, "C1'")
            pos_c4_input = get_coords(residue, "C4'")

            pos_c3 = np.zeros((combinations, xyz_coords))
            pos_c5 = np.zeros((combinations, xyz_coords))
            pos_c2 = np.zeros((combinations, xyz_coords))
            pos_o4 = np.zeros((combinations, xyz_coords))
            pos_o2 = np.zeros((combinations, xyz_coords))

            for i,variants in enumerate(all_variants):
                pos_c5[i], pos_o4[i], pos_c3[i], pos_c2[i], pos_o2[i] = infer_posits(
                    [pos_o5_input, pos_c1_input, pos_c4_input], variants
                )


            pos_o3 = np.repeat(pos_o3_input[np.newaxis,:], combinations, axis = 0)
            pos_c1 = np.repeat(pos_c1_input[np.newaxis,:], combinations, axis = 0)
            pos_c4 = np.repeat(pos_c4_input[np.newaxis,:], combinations, axis = 0)

            dist_c1_c2 = np.linalg.norm(pos_c1 - pos_c2, axis = 1)
            dist_c2_c3 = np.linalg.norm(pos_c2 - pos_c3, axis = 1)
            dist_c3_c4 = np.linalg.norm(pos_c3 - pos_c4, axis = 1)
            dist_c4_o4 = np.linalg.norm(pos_c4 - pos_o4, axis = 1)
            dist_o4_c1 = np.linalg.norm(pos_o4 - pos_c1, axis = 1)
            dist_c3_o3 = np.linalg.norm(pos_c3 - pos_o3, axis = 1)
            dist_o2_o3 = np.linalg.norm(pos_o2 - pos_o3, axis = 1)
            dist_o2_c3 = np.linalg.norm(pos_o2 - pos_c3, axis = 1)

            penalty = sum((
                np.abs(dist_c1_c2 - avg_c1_c2_dist),
                np.abs(dist_c2_c3 - avg_c2_c3_dist),
                np.abs(dist_c3_c4 - avg_c3_c4_dist),
                np.abs(dist_c4_o4 - avg_c4_o4_dist),
                np.abs(dist_o4_c1 - avg_o4_c1_dist),
                np.abs(dist_c3_o3 - avg_c3_o3_dist),
                np.abs(dist_o2_o3 - avg_o2_o3_dist),
                np.abs(dist_o2_c3 - avg_o2_c3_dist),
            ))
            idx_best = np.argmin(penalty)

            self._reconstruct_fa_particle(residue, "C3'", pos_c3[idx_best])
            self._reconstruct_fa_particle(residue, "C5'", pos_c5[idx_best])
            self._reconstruct_fa_particle(residue, "C2'", pos_c2[idx_best])
            self._reconstruct_fa_particle(residue, "O4'", pos_o4[idx_best])
            self._reconstruct_fa_particle(residue, "O2'", pos_o2[idx_best])


    # --------------------------------------------------------------------------
    def _reconstruct_nitro_base(self) -> None:
        nc: callable = h2f.Mapping.name_centroid

        def rec(residue: h2f.Residue, name_target: str, names_ref: tuple[str, str, str]):
            geo_ref = residue.init_fa_geometries(*names_ref)
            pos_new = self._get_new_pos(geo_ref, name_target, 0, restype = residue.restype)
            self._reconstruct_fa_particle(residue, name_target, pos_new)

        def reconstruct_u():
            rec(residue, "C6", (  "P", "C1'", nc(1)))
            rec(residue, "N1", ("C1'", nc(1),  "C6"))
            rec(residue, "C2", ( "C6",  "N1", nc(1)))
            rec(residue, "O2", ( "C6",  "N1", nc(1)))
            rec(residue, "N3", ( "C6",  "N1", nc(1)))
            rec(residue, "C4", ( "N1",  "C6", nc(1)))
            rec(residue, "O4", (nc(1),  "N3",  "C4"))
            rec(residue, "C5", ( "N1",  "C6", nc(1)))

        def reconstruct_c():
            rec(residue, "C6", ("C4'", "C1'", nc(1)))
            rec(residue, "N1", ("C1'", nc(1),  "C6"))
            rec(residue, "C2", ( "C6",  "N1", nc(1)))
            rec(residue, "O2", ( "C6",  "N1", nc(1)))
            rec(residue, "N3", ( "C6",  "N1", nc(1)))
            rec(residue, "C4", ( "N1",  "C6", nc(1)))
            rec(residue, "N4", (nc(1),  "N3",  "C4"))
            rec(residue, "C5", ( "N1",  "C6", nc(1)))

        def reconstruct_a():
            rec(residue, "N9", ("C1'", nc(2), nc(1)))
            rec(residue, "C4", ( "N9", nc(2), nc(1)))
            rec(residue, "N3", ( "N9", nc(1), nc(2)))
            rec(residue, "C2", ( "N9", nc(1), nc(2)))
            rec(residue, "N1", ( "N9", nc(1), nc(2)))
            rec(residue, "C6", ( "N9", nc(1), nc(2)))
            rec(residue, "N6", ( "N9", nc(1), nc(2)))
            rec(residue, "C5", ( "N9", nc(1), nc(2)))
            rec(residue, "N7", ( "N9", nc(2), nc(1)))
            rec(residue, "C8", (nc(2), nc(1),  "N9"))

        def reconstruct_g():
            rec(residue, "N9", ("C1'", nc(2), nc(1)))
            rec(residue, "C4", ( "N9", nc(2), nc(1)))
            rec(residue, "N3", ( "N9", nc(1), nc(2)))
            rec(residue, "C2", ( "N9", nc(1), nc(2)))
            rec(residue, "N2", ( "N9", nc(1), nc(2)))
            rec(residue, "N1", ( "N9", nc(1), nc(2)))
            rec(residue, "C6", ( "N9", nc(1), nc(2)))
            rec(residue, "O6", ( "N9", nc(1), nc(2)))
            rec(residue, "C5", ( "N9", nc(1), nc(2)))
            rec(residue, "N7", ( "N9", nc(2), nc(1)))
            rec(residue, "C8", (nc(2), nc(1),  "N9"))

        for residue in self.residues:
            match residue.restype:
                case 'U': reconstruct_u()
                case 'C': reconstruct_c()
                case 'A': reconstruct_a()
                case 'G': reconstruct_g()


    # --------------------------------------------------------------------------
    def _drop_temp_particles(self) -> None:
        for residue in self.residues:
            residue.fa_data.pop_filtered(atomname = h2f.Mapping.name_centroid(1))
            if residue.restype not in ("A", "G"): continue
            residue.fa_data.pop_filtered(atomname = h2f.Mapping.name_centroid(2))


    # --------------------------------------------------------------------------
    def _get_new_pos(self,
        geo_ref: h2f.Geometry, fa_name_target: str,
        idx_variant: int, restype: str = '_',
    ) -> np.ndarray:

        this_model = self.model[restype][fa_name_target]
        if geo_ref is None: return

        geo_ref.set_geo(*this_model[idx_variant])
        return geo_ref.infer_fourth_particle()


    # --------------------------------------------------------------------------
    def _reconstruct_fa_particle(self,
        residue_orig: h2f.Residue,
        fa_name_target: str,
        new_pos: np.ndarray | None,
    ) -> None:
        if new_pos is None: return
        particle = h2f.PDBTable(sections = {
            "atomid"  : (str(len(residue_orig.fa_data)),),
            "atomname": (fa_name_target,),
            "element" : (fa_name_target[0],),
            "resname" : (residue_orig.resname,),
            "chainid" : (residue_orig.chainid,),
            "segid"   : (residue_orig.chainid,),
            "resid"   : (residue_orig.resid,),
            "xcoords" : (f"{new_pos[0]:.3f}",),
            "ycoords" : (f"{new_pos[1]:.3f}",),
            "zcoords" : (f"{new_pos[2]:.3f}",),
        })
        residue_orig.fa_data.append_table(particle)


# //////////////////////////////////////////////////////////////////////////////
