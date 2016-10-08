import os
import argparse


if __name__ == '__main__':

    parser = argparse.ArgumentParser()
    parser.add_argument('--xtc', help='trajectory filename')
    parser.add_argument('--system', help='system name')
    parser.add_argument('--xtc_nopbc', help='trajectory without PBC')
    parser.add_argument('--reference', help='reference pdb name')
    parser.add_argument('--lh5', help='lh5 filename')
    parser.add_argument('--tpr', help='tpr filename')
    args = parser.parse_args()

    ## Catch output of gromacs trajcetory
    os.system('gmx check -f {0} 2> error_stream.log'.format(args.xtc))


    ## Feed output + check the trajectory, get the traj info
    os.system('cat error_stream.log | python checktrajectory.py')

    # Remove SOL and PBC
    os.system('echo {0} | gmx trjconv -f {1} -s {3} -o {2} -pbc mol'.format(args.system, args.xtc, args.xtc_nopbc, args.tpr))

    # Convert xtc to lh5    
    os.system('python convert2lh5.py --inp {0} --ref {1} --out {2}'.format(args.xtc_nopbc, args.reference, args.lh5))