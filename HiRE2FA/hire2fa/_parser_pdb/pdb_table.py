import numpy as np
from pathlib import Path

import hire2fa as h2f

# //////////////////////////////////////////////////////////////////////////////
class PDBTable:
    def __init__(self, sections: dict[str, tuple[str]] = None):
        self.sections = {
            name: tuple() for name in h2f.PDBConstants.SECTION_NAMES
        } if sections is None else sections.copy()

    # --------------------------------------------------------------------------
    def __len__(self): return len(self.sections["atomid"])

    # --------------------------------------------------------------------------
    def __repr__(self):
        sorted_keys = ["atomid", "atomname", "resname", "chainid", "resid", "xcoords", "ycoords", "zcoords"]
        lines = zip(*(self.sections[k] for k in sorted_keys))
        header = '\n'.join((
            f"PDBTable with {len(self)} particles",
            ''.join(x.ljust(8) for x in sorted_keys)
        )) + '\n'
        body = '\n'.join(
            ''.join(x.ljust(8) for x in line) for line in lines
        ) + '\n'
        return header + body


   # --------------------------------------------------------------------------
    @classmethod
    def join(cls, tables: list["PDBTable"]) -> "PDBTable":
        """Joins multiple `PDBTable` objects into a single one. The order of the sections is preserved."""
        out = cls()
        for table in tables:
            out.append_table(table)
        return out


   # --------------------------------------------------------------------------
    @classmethod
    def read_pdb(cls, path_pdb: str | Path):
        data = list(filter(
            lambda line: line.startswith("ATOM"),
            Path(path_pdb).read_text().strip().splitlines()
        ))
        ### current sections: atomid / atomname / resname / chainid / resid / xcoords / ycoords / zcoords
        ### when exporting...:atomid / atomname / resname / chainid / resid / xcoords / ycoords / zcoords / segid / element
        sections: dict[str, tuple[str]] = {
            name: tuple(h2f.PDBConstants.extract_section(line, name) for line in data)
            for name in h2f.PDBConstants.SECTION_NAMES
        }
        return cls(sections)

    # --------------------------------------------------------------------------
    def write_pdb(self, path_pdb: str | Path) -> None:
        out = '\n'.join(
            h2f.PDBConstants.export_row(self.sections, i) for i in range(len(self))
        )
        with open(path_pdb, 'w') as file:
            file.write(out)


    # --------------------------------------------------------------------------
    def append_line(self, line: str) -> None:
        """Appends an extra element for every section. Data is extracted from the `line` provided."""
        for section in self.sections.keys():
            self.sections[section] += (h2f.PDBConstants.extract_section(line, section),)

    # --------------------------------------------------------------------------
    def append_table(self, other: "PDBTable") -> None:
        """Appends another `PDBTable` to this one."""
        for section,other_tup in other.sections.items():
            self.sections[section] += other_tup


    # --------------------------------------------------------------------------
    def copy(self) -> "PDBTable":
        """Returns a copy of this `PDBTable`."""
        return PDBTable(sections = self.sections)

    # --------------------------------------------------------------------------
    def copy_filtered(self, **filter_kvs) -> "PDBTable":
        """Returns an unordered copy of this `PDBTable` that matches the filter values for the specified sections"""
        idxs = set(range(len(self)))
        for section,val in filter_kvs.items():
            idxs.intersection_update(self._get_filter_idxs(section, val))

        copy = PDBTable(sections = self.sections.copy())
        copy.slice_sections(tuple(idxs))
        return copy

    # --------------------------------------------------------------------------
    def pop_filtered(self, **filter_kvs) -> "PDBTable":
        """
        Returns an unordered copy of this `PDBTable` that matches the filter values for the specified sections.
        Removes the matched rows from the original `PDBTable` instance.
        """
        idxs_all = set(range(len(self)))
        idxs = idxs_all.copy()
        for section,val in filter_kvs.items():
            idxs.intersection_update(self._get_filter_idxs(section, val))

        copy = PDBTable(sections = self.sections.copy())
        copy.slice_sections(tuple(idxs))
        self.slice_sections(tuple(idxs_all - idxs))
        return copy


    # --------------------------------------------------------------------------
    def normalize_indices(self) -> None:
        """Normalizes the `atomid` section to be a sequence of integers starting from 1. This is useful after filtering or joining `PDBTable` objects."""
        self.sections["atomid"] = tuple(str(i) for i in range(1, len(self) + 1))


    # --------------------------------------------------------------------------
    def sort(self) -> None:
        """Sort every section by using the `atomid` section as sorting reference"""
        self.slice_sections(idxs = h2f.Utils.argsort(self.sections["atomid"]))


    # --------------------------------------------------------------------------
    def slice_sections(self, idxs: list[int] | tuple[int]):
        """In-place slice of every section by using a list of indices"""
        for k,tup in self.sections.items():
            self.sections[k] = tuple(tup[i] for i in idxs)


    # --------------------------------------------------------------------------
    def safe_get_first_value(self, section: str) -> str | None:
        """Returns the first value of a section if it exists, otherwise returns `None`"""
        return self.sections[section][0] if len(self) else None

    # --------------------------------------------------------------------------
    def get_coords_array(self) -> np.ndarray:
        """Returns the coordinates of this `PDBTable` as a numpy array of shape (N, 3)"""
        return np.array(list(zip(
            map(float, self.sections["xcoords"]),
            map(float, self.sections["ycoords"]),
            map(float, self.sections["zcoords"]),
        )))


    # --------------------------------------------------------------------------
    def iter_chains(self):
        """Yields `PDBTable` objects for each chain in this `PDBTable`. Chains are identified and **sorted** by their `chainid` section."""
        for chainid in sorted(set(self.sections["chainid"])):
            yield PDBTable.copy_filtered(self, chainid = chainid)

    # --------------------------------------------------------------------------
    def iter_residues(self):
        """Yields `PDBTable` objects for each residue in this `PDBTable`. Residues are identified and **sorted** by their `resid` section."""
        for resid in sorted(set(self.sections["resid"]), key = int):
            yield PDBTable.copy_filtered(self, resid = resid)


    # --------------------------------------------------------------------------
    def drop_duplicates(self) -> None:
        """In-place drops duplicate rows from this `PDBTable` by using the `atomname` section as reference. Useful for removing unwanted altlocs"""
        seen = set()
        idxs = []
        for i,atomname in enumerate(self.sections["atomname"]):
            if atomname in seen: continue
            seen.add(atomname)
            idxs.append(i)
        self.slice_sections(idxs)


    # --------------------------------------------------------------------------
    def _get_filter_idxs(self, section: str, value: str) -> list[int]:
        """Returns the row indices for atoms that match a given `value` for a PDB `section`"""
        return [i for i,v in enumerate(self.sections[section]) if v == value]


# //////////////////////////////////////////////////////////////////////////////
