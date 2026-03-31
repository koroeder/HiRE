import hire2fa as h2f

# //////////////////////////////////////////////////////////////////////////////
class Residue:
    def __init__(self, data: h2f.PDBTable, init_with_cg_data: bool = True):
        assert len(data) > 0, "Empty residue"

        self.resname: str = data.safe_get_first_value("resname").upper()
        self.resid  : str = data.safe_get_first_value("resid")
        self.chainid: str = data.safe_get_first_value("chainid")

        data = self._assert_generic_residue(data)

        self.is_5terminus = self.resname.endswith("5")
        self.restype: str = self.resname[0].upper() # U / C / A / G

        self.cg_bead_names = h2f.Mapping.EXPECTED_CG_BEADS[self.restype].copy()
        if self.is_5terminus: self.cg_bead_names.remove("P")

        if init_with_cg_data:
            self.fa_data = h2f.PDBTable()
            self.cg_data = data
            self._assert_cg_residue()
        else:
            self.fa_data = data
            self.cg_data = h2f.PDBTable()

        self.fa_data.sort()
        self.cg_data.sort()


    # --------------------------------------------------------------------------
    @classmethod
    def get_pdb_residues(cls, pdb_table: h2f.PDBTable, init_with_cg_data: bool = True) -> list["Residue"]:
        """
        Splits a `h2f.PDBTable` into a list of `Residue` objects, one for each chain and residue.
        The residues are sorted by their `chainid` and `resid` sections.
        """
        return [
            Residue(residue, init_with_cg_data)
            for chain in pdb_table.iter_chains()
            for residue in chain.iter_residues()
        ]

    # --------------------------------------------------------------------------
    @classmethod
    def get_pdb_chains(cls, pdb_table: h2f.PDBTable, init_with_cg_data: bool = True) -> list[list["Residue"]]:
        """
        Splits a `h2f.PDBTable` into a list of lists of `Residue` objects.
        A list of residues represents a different chain.
        The residues are sorted by their `chainid` and `resid` sections.
        """
        return [
            [
                Residue(residue, init_with_cg_data)
                for residue in chain.iter_residues()
            ]
            for chain in pdb_table.iter_chains()
        ]


    # --------------------------------------------------------------------------
    @classmethod
    def apply_o3_shift(cls, chains: list[list["Residue"]], do_cg: bool) -> None:
        """Moves O3' atoms from the end of each residue to the beginning of the next one."""
        ##### PART 0: make sure indices are normalized
        for chain in chains:
            for residue in chain:
                data = residue.get_data(do_cg)
                data.normalize_indices()

        ##### PART 1: Move O3' atoms between neighbouring residues
        name_o3 = "O3" if do_cg else "O3'"

        for chain in chains:
            if len(chain) < 2: continue

            for i,(res_this,res_prev) in enumerate(zip(chain[:0:-1], chain[-2::-1])):
                res_this: Residue
                res_prev: Residue
                data_this = res_this.get_data(do_cg)
                data_prev = res_prev.get_data(do_cg)

                if not i: # no need for the O3 atom of the chain's last residue
                    data_this.pop_filtered(atomname = name_o3)

                atom_o3 = data_prev.pop_filtered(atomname = name_o3)

                if not len(atom_o3): continue
                if len(atom_o3) > 1:
                    raise ValueError(f"Multiple atoms with the name {name_o3} in the same residue")

                resid = data_this.safe_get_first_value("resid")
                atom_o3.sections["resid"] = (resid,)

                data_this.append_table(atom_o3)


    # --------------------------------------------------------------------------
    def get_data(self, do_cg: bool) -> h2f.PDBTable:
        """Returns either `cg_data` or `fa_data` depending on the value of `do_cg`."""
        return self.cg_data if do_cg else self.fa_data

    # --------------------------------------------------------------------------
    def get_functional_group(self, atomnames: list[str]) -> h2f.PDBTable:
        return h2f.PDBTable.join(
            self.fa_data.copy_filtered(atomname = name)
            for name in atomnames
        )


    # --------------------------------------------------------------------------
    def copy(self) -> "Residue":
        """Returns a copy of this `Residue`."""
        res_copy = Residue.__new__(Residue)
        res_copy.resname = self.resname
        res_copy.resid = self.resid
        res_copy.chainid = self.chainid
        res_copy.restype = self.restype
        res_copy.is_5terminus = self.is_5terminus
        res_copy.cg_bead_names = self.cg_bead_names.copy()
        res_copy.fa_data = self.fa_data.copy()
        res_copy.cg_data = self.cg_data.copy()
        return res_copy


    # --------------------------------------------------------------------------
    def init_fa_geometries(self, name_0: str, name_1: str, name_2: str) -> h2f.Geometry:
        tab_0 = self.fa_data.copy_filtered(atomname = name_0)
        tab_1 = self.fa_data.copy_filtered(atomname = name_1)
        tab_2 = self.fa_data.copy_filtered(atomname = name_2)
        if not (len(tab_0) == len(tab_1) == len(tab_2) == 1): return

        coords_0 = tab_0.get_coords_array()[0]
        coords_1 = tab_1.get_coords_array()[0]
        coords_2 = tab_2.get_coords_array()[0]
        return h2f.Geometry(coords_0, coords_1, coords_2)


    # --------------------------------------------------------------------------
    def calc_geometry_values(self, geometry: h2f.Geometry | None, name: str) -> h2f.Geometry | None:
        if geometry is None: return

        tab = self.fa_data.copy_filtered(atomname = name)
        if len(tab) != 1: return # atom name isn't unique

        coords = tab.get_coords_array()[0]
        return geometry.calc_geo(coords)


    # --------------------------------------------------------------------------
    def _assert_generic_residue(self, data: h2f.PDBTable) -> h2f.PDBTable:
        if any(r != self.resname for r in data.sections["resname"]):
            raise ValueError("Inconsistent residue names in the same residue")

        if self.resname not in h2f.Mapping.KNOWN_RESIDUES:
            raise ValueError(f"Unknown residue name: {self.resname}")

        particle_names = set(data.sections["atomname"])
        if len(data.sections["atomname"]) > len(particle_names):
            print("WARNING: duplicate particle names in the same residue. Dropping the duplicates.")
            data.drop_duplicates()

        return data


    # --------------------------------------------------------------------------
    def _assert_cg_residue(self):
        bead_names = set(self.cg_data.sections["atomname"])
        expected_beads = self.cg_bead_names
        if bead_names != expected_beads:
            raise ValueError(f"Unexpected bead names for residue {self.resname} (found: {bead_names}, expected: {expected_beads})")


# //////////////////////////////////////////////////////////////////////////////
