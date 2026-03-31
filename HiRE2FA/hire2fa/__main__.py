import sys
from pathlib import Path

import hire2fa as h2f

# ------------------------------------------------------------------------------
def print_usage_and_exit(exit_code: int):
    print(f"HiRE2FA (v{h2f.__version__}). Usage:")
    print(f"    {sys.argv[0]} [path_in_cg] [path_out_fa]")
    exit(exit_code)

# ------------------------------------------------------------------------------
def main():
    if len(sys.argv) < 3: print_usage_and_exit(-1)
    if sys.argv[1] in ("-h", "--help"): print_usage_and_exit(0)
    if sys.argv[2] in ("-h", "--help"): print_usage_and_exit(0)

    PATH_PDB_CG = Path(sys.argv[1])
    PATH_PDB_FA = Path(sys.argv[2])

    rec = h2f.Reconstructor(
        path_pdb_cg = PATH_PDB_CG,
        path_model = Path(__file__).parent / "_data/model.json"
    )
    rec.reconstruct()
    rec.export_reconstructed(PATH_PDB_FA)


################################################################################
if __name__ == "__main__":
    main()


################################################################################
