import mdtraj
import sys
import pylab
import matplotlib as mpl
import matplotlib.pyplot as plt

if __name__ == '__main__':

    fname = sys.argv[1]

    traj = mdtraj.load(fname)
    atoms, bonds = traj.topology.to_dataframe()
    print atoms

    psi_indices, phi_indices = [6, 8, 14, 16], [4, 6, 8, 14]
    angles = mdtraj.geometry.compute_dihedrals(traj, [phi_indices, psi_indices])

    plt.scatter(angles[:, 0], angles[:, 1], marker='x', c=traj.time)
