import argparse
from os.path import exists

def create_start(coords):
    outf = open("start","w")
    for xyz in coords:
        line = "{:16.10f}".format(xyz[0]) + "{:16.10f}".format(xyz[1]) + "{:16.10f}".format(xyz[2]) + "\n"
        outf.write(line)
    outf.close()

def parse_pdb(pdb_file):
    coords = list()
    with open(pdb_file, "r") as f:
        for line in f:
            if line[0:4]=="ATOM":
                xyz = line[30:54]
                xyz = [float(i) for i in xyz.split()]
                coords.append((xyz[0],xyz[1],xyz[2]))
    return coords

if __name__ == "__main__":
    parser = argparse.ArgumentParser(prog = 'Start file creator',
                    description = 'This script creates a start file for HiRE from a CG pdb file',
                    epilog = 'Do not forget that you also need a topology file. If you create a topology from the pdb, a start file will be created automatically.')

    parser.add_argument('filename', help='Name of CG pdb file used to create start file') 
    args = parser.parse_args()
    fname = args.filename
    if not(exists(fname)):
        raise FileNotFoundError("File %s does not exist" % fname)
    create_start(parse_pdb(fname))
    