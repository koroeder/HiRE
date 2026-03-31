# //////////////////////////////////////////////////////////////////////////////
class Mapping:
    KNOWN_RESIDUES = {
        "U", "U3", "U5",
        "C", "C3", "C5",
        "A", "A3", "A5",
        "G", "G3", "G5",
    }
    EXPECTED_CG_BEADS = {
        'U': {"P", "O3", "O5", "R1", "R4", "U1"},
        'C': {"P", "O3", "O5", "R1", "R4", "C1"},
        'A': {"P", "O3", "O5", "R1", "R4", "A1", "A2"},
        'G': {"P", "O3", "O5", "R1", "R4", "G1", "G2"},
    }
    ATOM_MASSES = {
        "C": 12.011, "H": 1.008, "O": 16.000, "N": 14.007, "P": 30.974,
    }
    BASENAMES_FOR_CENTROID_1 = {
        'U': ["N1", "C2", "O2", "N3", "C4", "O4", "C5", "C6"],
        'C': ["N1", "C2", "O2", "N3", "C4", "N4", "C5", "C6"],
        'A': ["C4", "C5", "N7", "C8", "N9"],
        'G': ["C4", "C5", "N7", "C8", "N9"],
    }
    BASENAMES_FOR_CENTROID_2 = {
        'A': ["C4", "C5", "C6", "N1", "C2", "N3"],
        'G': ["C4", "C5", "C6", "N1", "C2", "N3"],
    }
    MAP_STATSNAME_TO_STANDARDNAME = {
        "op1": "OP1", "op2": "OP2",
        "c5": "C5'", "c3": "C3'", "c2": "C2'", "o2": "O2'", "o4": "O4'",
        "buC6" : "C6", "buN1" : "N1", "buC2" : "C2", "buO2" : "O2", "buN3" : "N3", "buC4" : "C4", "buO4" : "O4", "buC5" : "C5",
        "bcC6" : "C6", "bcN1" : "N1", "bcC2" : "C2", "bcO2" : "O2", "bcN3" : "N3", "bcC4" : "C4", "bcN4" : "N4", "bcC5" : "C5",
        "baN9" : "N9", "baC4" : "C4", "baN3" : "N3", "baC2" : "C2", "baN1" : "N1", "baC6" : "C6", "baN6" : "N6", "baC5" : "C5", "baN7" : "N7", "baC8" : "C8",
        "bgN9" : "N9", "bgC4" : "C4", "bgN3" : "N3", "bgC2" : "C2", "bgN2" : "N2", "bgN1" : "N1", "bgC6" : "C6", "bgO6" : "O6", "bgC5" : "C5", "bgN7" : "N7", "bgC8" : "C8",
    }
    _DIRECT_MAPPING_GENERIC = {
        ( "P",   "P"),
        ("O3", "O3'"),
        ("O5", "O5'"),
        ("R1", "C1'"),
        ("R4", "C4'"),
    }

    # --------------------------------------------------------------------------
    def name_centroid(one_or_two: int | str):
        return f"_BC{one_or_two}"

    # --------------------------------------------------------------------------
    @classmethod
    def iter_direct_mapping(cls):
        yield from cls._DIRECT_MAPPING_GENERIC
        yield "U1", cls.name_centroid(1)
        yield "C1", cls.name_centroid(1)
        yield "A1", cls.name_centroid(1)
        yield "G1", cls.name_centroid(1)
        yield "A2", cls.name_centroid(2)
        yield "G2", cls.name_centroid(2)


# //////////////////////////////////////////////////////////////////////////////
