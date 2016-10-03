import argparse
import msmbuilder.Trajectory

def convertXtc2lh5(xtcfile,  ref_conf, outFilename):
    ''' Convert the xtc-files to a .lh5 file '''

    Traj = msmbuilder.Trajectory.Trajectory.LoadFromXTC([xtcfile], PDBFilename=ref_conf)
    Traj.Save("%s"%outFilename)

if __name__ == '__main__':

    parser = argparse.ArgumentParser()
    parser.add_argument("--inp", help="input filename")
    parser.add_argument("--out", help="output filename")
    parser.add_argument("--ref", help="reference filename")
    args = parser.parse_args()

    convertXtc2lh5(args.inp, args.ref, args.out)