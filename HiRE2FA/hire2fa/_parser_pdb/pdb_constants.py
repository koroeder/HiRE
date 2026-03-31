# //////////////////////////////////////////////////////////////////////////////
class PDBConstants:
    _PDB_IDXS: dict[str, int] = dict(
        atomid_START = 6,
        atomid_END   = 11,
        atomname_START = 12,
        atomname_END   = 16,
        resname_START = 17,
        resname_END   = 20,
        chainid_START = 21,
        chainid_END   = 22,
        resid_START = 22,
        resid_END   = 26,
        xcoords_START = 30,
        xcoords_END   = 38,
        ycoords_START = 38,
        ycoords_END   = 46,
        zcoords_START = 46,
        zcoords_END   = 54,
        segid_START = 72,
        segid_END   = 76,
        element_START = 76,
        element_END   = 78,
    )
    _PDB_JUSTIFICATIONS: dict[str, callable] = dict(
        atomid   = str.rjust,
        atomname = str.ljust,
        resname  = str.rjust,
        chainid  = str.ljust,
        resid    = str.rjust,
        xcoords  = str.rjust,
        ycoords  = str.rjust,
        zcoords  = str.rjust,
        segid    = str.ljust,
        element  = str.rjust,
    )
    _LENGTH_RECORD: int = 80
    SECTION_NAMES = set(_PDB_JUSTIFICATIONS.keys())

    # --------------------------------------------------------------------------
    @classmethod
    def extract_section(cls, line: str, name: str) -> str:
        start = cls._PDB_IDXS[f"{name}_START"]
        end   = cls._PDB_IDXS[f"{name}_END"]
        return line[start:end].strip()


    # --------------------------------------------------------------------------
    @classmethod
    def export_row(cls, data: dict[str, tuple[str]], row_idx: int) -> str:
        chars: list[str] = [' ' for _ in range(cls._LENGTH_RECORD)]
        chars[:6] = "ATOM  "

        for name in cls.SECTION_NAMES:
            start: int = cls._PDB_IDXS[f"{name}_START"]
            end  : int = cls._PDB_IDXS[f"{name}_END"]
            just : callable  = cls._PDB_JUSTIFICATIONS[name]
            width: int = end - start
            value: str = just(data[name][row_idx], width)[:width]
            chars[start:end] = value

        return ''.join(chars)


    # --------------------------------------------------------------------------
# //////////////////////////////////////////////////////////////////////////////
