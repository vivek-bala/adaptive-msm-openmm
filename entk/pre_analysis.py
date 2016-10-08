import os


if __name__ == '__main__':


	## Catch output of gromacs trajcetory
	os.system('gmx check -f *.xtc 2> error_stream.log')


	## Feed output + check the trajectory, get the traj info
	os.system('cat error_stream.log | python checktrajectory.py')

	# Remove SOL and PBC
    os.system('echo $system | gmx trjconv -f *.xtc -s *.tpr -o traj_$inst.nopbc.xtc -pbc mol')