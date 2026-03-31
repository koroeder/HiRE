import numpy as np

import hire2fa as h2f

# //////////////////////////////////////////////////////////////////////////////
class Geometry:
    def __init__(self, pos_a: np.ndarray, pos_b: np.ndarray, pos_c: np.ndarray):
        self.pos_a = pos_a.copy() # aka "a"
        self.pos_b = pos_b.copy() # aka "b"
        self.pos_c = pos_c.copy() # aka "c"
        self.dist : float | None = None # aka "length" or "cd"
        self.angle: float | None = None # aka "theta" or "bcd"
        self.dihed: float | None = None # aka "phi" or "abcd"


    # --------------------------------------------------------------------------
    @staticmethod
    def calc_center_of_mass(coords: np.ndarray, names = list[str]) -> np.ndarray:
        masses = np.array([h2f.Mapping.ATOM_MASSES[name[0]] for name in names])
        return np.sum(
            coords * masses[:, np.newaxis],
            axis = 0
        ) / np.sum(masses)


    # --------------------------------------------------------------------------
    def copy(self) -> "Geometry":
        copy = Geometry(self.pos_a, self.pos_b, self.pos_c)
        copy.dist = self.dist
        copy.angle = self.angle
        copy.dihed = self.dihed
        return copy


    # --------------------------------------------------------------------------
    def calc_geo(self, pos_d: np.ndarray) -> "Geometry":
        out = self.copy()
        out._calc_dist (pos_d)
        out._calc_angle(pos_d)
        out._calc_dihed(pos_d)
        return out


    # --------------------------------------------------------------------------
    def set_geo(self, dist: float, angle: float, dihed: float) -> None:
        self.dist = dist
        self.angle = angle
        self.dihed = dihed


    # --------------------------------------------------------------------------
    def infer_fourth_particle(self) -> np.ndarray:
        """
        Returnns a fourth position `d` in cartesian coordinates, by using as reference the three positions `a`, `b`, `c`
        and the geometry values dist `cd`, angle `bcd` and dihedral `abcd`.
        The geometry values need to be set beforehand by calling `set_geo` or `calc_geo`.
        Implemention of SN-NeRF, described in:

        > Parsons, J., Holmes, J. B., Rojas, J. M., Tsai, J., & Strauss, C. E. M. (2005).
        Practical conversion from torsion space to Cartesian space for in silico protein synthesis.
        Journal of computational chemistry, 26(10), 1063-8. doi: 10.1002/jcc.20237.

        """

        if self.dist is None:
            raise ValueError("either `calc_geo` or `set_geo` needs to be called before infer_fourth_particle")

        a = self.pos_a
        b = self.pos_b
        c = self.pos_c

        d0 = self.dist * np.array([
            np.cos(self.angle),
            np.sin(self.angle) * np.cos(self.dihed),
            np.sin(self.angle) * np.sin(self.dihed),
        ])

        ab = b - a
        bc = c - b
        bc /= np.linalg.norm(bc)

        n = np.cross(ab, bc)
        n /= np.linalg.norm(n)

        nxbc = np.cross(n, bc)
        M = np.array([bc, nxbc, n]).T
        dstar = M @ d0
        return c + dstar


    # --------------------------------------------------------------------------
    def _calc_dist(self, pos_d: np.ndarray) -> None:
        self.dist = np.linalg.norm(pos_d - self.pos_c)


    # --------------------------------------------------------------------------
    def _calc_angle(self, pos_d: np.ndarray) -> None:
        # Assumes the following order of atoms. Note that the angle is the one between the BC axis and the point D
        # It's NOT the angle between the points BCD.
        #        D
        #       / (angle)
        # B -- C -- axis_bc

        axis_bc = self.pos_c - self.pos_b
        vec_cd = pos_d - self.pos_c
        cosine = np.dot(axis_bc, vec_cd) / (np.sqrt(np.dot(axis_bc, axis_bc) * np.dot(vec_cd, vec_cd)))
        self.angle = np.arccos(cosine) # in radians


    # --------------------------------------------------------------------------
    def _calc_dihed(self, pos_d: np.ndarray) -> None:
        # Assumes the following order of atoms for proper dihedrals:
        #   B   D :   B B   D
        #  / \ /  :  / \ \ /
        # A   C   : A   C C

        vec_ab = self.pos_b - self.pos_a
        vec_bc = self.pos_c - self.pos_b
        vec_cd = pos_d - self.pos_c

        n1 = np.cross(vec_ab, vec_bc)
        n2 = np.cross(vec_bc, vec_cd)
        cosine = np.dot(n1, n2) / (np.sqrt(np.dot(n1, n1) * np.dot(n2, n2)))

        sign = 1 if np.dot(vec_bc, np.cross(n1, n2)) > 0 else -1
        self.dihed = sign * np.arccos(np.clip(cosine, -1.0, 1.0)) # in radians


# //////////////////////////////////////////////////////////////////////////////
