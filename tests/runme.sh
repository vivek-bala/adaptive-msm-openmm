#!/bin/sh

## Setup GROMACS
export PATH=$PATH:/home/vivek/Research/modules/gromacs-5.1.3/build/bin
export PYTHONPATH=$PYTHONPATH:/home/vivek/Research/repos/adaptive-msm/tests

## Source location, cur location
src=/home/vivek/Research/repos/adaptive-msm
cur=$PWD

## Input parameters
system=Protein
sims=4

## Create shared space
rm shared -rf
mkdir shared
shared=$cur/shared

#-----------------------------------------------------------------------------------------
## Preprocessing stage
#-----------------------------------------------------------------------------------------
# Create folder for first stage
cd $cur
rm preproc -rf
mkdir preproc
cd preproc
preproc=$cur/preproc

# Copy configuration and topology files
cp $src/*.gro .
cp $src/*.top . 
cp $src/*.mdp .

# Run preprocessing kernel on all configuration files
inst=0
for file in ./*.gro
do 
    mkdir sim_$inst
    cd sim_$inst
    cp ../equil$inst.gro .
    cp ../topol.top .
    cp ../grompp.mdp .
    gmx grompp -f grompp.mdp -quiet -c $file -p topol.top -o topol.tpr
    cp topol.tpr $shared/topol_$inst.tpr
    cd ..
    inst=$((inst+1))
done
#-----------------------------------------------------------------------------------------


#-----------------------------------------------------------------------------------------
## Generate pdb files
#-----------------------------------------------------------------------------------------
# Create folder for second stage
cd $cur
rm gen_pdb -rf
mkdir gen_pdb
cd gen_pdb

# Copy files to current folder
cp $shared/*.tpr .
cp $src/*.gro .
cp $src/generate_reference.py .

# Generate pdbs
inst=0
for file in ./*.tpr
do
    mkdir sim_$inst
    cd sim_$inst
    cp ../topol_$inst.tpr .
    cp ../equil$inst.gro .
    cp ../generate_reference.py .
    echo $system | gmx trjconv -f equil$inst.gro -s topol_$inst.tpr -o ref_tmp.pdb
    python generate_reference.py --inp ref_tmp.pdb --out reference.pdb
    cp reference.pdb $shared/reference_$inst.pdb
    cd ..
    inst=$((inst+1))
done
#-----------------------------------------------------------------------------------------

#-----------------------------------------------------------------------------------------
## Tune and run initial simulations
#-----------------------------------------------------------------------------------------
# Create folder for third stage
cd $cur
rm tune_mdrun -rf
mkdir tune_mdrun
cd tune_mdrun

# Copy files to current folder
cp $src/*.gro .
cp $src/*.mdp .
cp $src/*.top .


# Perform grompp again
inst=0
for file in ./*.gro
do
    mkdir sim_$inst
    cd sim_$inst
    cp ../equil$inst.gro .
    cp ../topol.top .
    cp ../grompp.mdp .
    gmx grompp -f grompp.mdp -quiet -c $file -p topol.top -o topol.tpr
    cd ..
    inst=$((inst+1))
done


# Perform mdrun again
inst=0
for file in ./*.gro
do
    cd sim_$inst
    gmx mdrun -nt 1 -s topol.tpr -rcon 0.7 -maxh 0.0005
    cp mdout.mdp $shared/mdout_$inst.mdp
    cd ..
    inst=$((inst+1))
done
#-----------------------------------------------------------------------------------------

#-----------------------------------------------------------------------------------------
## Run simulation and perform trajectory collection

# This stage runs 
#-----------------------------------------------------------------------------------------
# Create folder for fourth stage
cd $cur
rm md_sims -rf
mkdir md_sims
cd md_sims

# Copy required files to current folder
cp $src/*.gro .
cp $src/*.top .
cp $src/procsettings.py .
cp $src/checktrajectory.py .
cp $src/convert2lh5.py .
cp $shared/mdout_*.mdp .
cp $shared/reference_*.pdb .

# Perform the simulation
for inst in `seq 0 $((sims-1))`;
do 

    # Create separate folder in shared 
    mkdir $shared/sim_$inst

    mkdir sim_$inst
    cd sim_$inst
    cp ../mdout_$inst.mdp .
    cp ../procsettings.py .
    
    # Create new mdp file
    python procsettings.py --inp mdout_$inst.mdp --out grompp.mdp 

    # Preprocessing 
    cp ../topol.top .
    cp ../equil$inst.gro .
    gmx grompp -f grompp.mdp -quiet -c equil$inst.gro -p topol.top -o topol.tpr


    # MD Run
    gmx mdrun -quiet -s topol.tpr -noappend -cpi state.cpt -rcon 0.7

    # Trajectory collection section
    gmx check -f *.xtc 2> error_stream.log
    cp ../checktrajectory.py .
    cat error_stream.log | python checktrajectory.py

    # Remove SOL and PBC
    echo $system | gmx trjconv -f *.xtc -s *.tpr -o traj_$inst.nopbc.xtc -pbc mol

    # Convert xtc to lh5
    cp ../reference_$inst.pdb .
    cp ../convert2lh5.py .
    python convert2lh5.py --inp traj_$inst.nopbc.xtc --ref reference_$inst.pdb --out file_$inst.lh5


    # Move back required data into shared
    cp file_$inst.lh5 $shared/sim_$inst/
    cp traj_$inst.nopbc.xtc $shared/sim_$inst/
    cp traj_comp.*.xtc $shared/sim_$inst/traj_$inst.xtc
    cp topol.tpr $shared/sim_$inst/topol_$inst.tpr
    cp traj_info.txt $shared/sim_$inst/traj_info_$inst.txt


    cd ..
done
#-----------------------------------------------------------------------------------------



#-----------------------------------------------------------------------------------------
## Run MSMproject
#-----------------------------------------------------------------------------------------
# Create folder for fifth stage
cd $cur
mkdir msm
cd msm

# Copy required files
cp $shared/sim_*/file_*.lh5 .
cp $shared/sim_*/traj_*.nopbc.xtc .
cp $shared/sim_*/traj_*.xtc .
cp $shared/sim_*/*.tpr .
cp $shared/sim_*/*.txt .
cp $shared/reference_0.pdb .
cp $src/MSMproject.py .

# Run the MSM kernel
python MSMproject.py --micro 100 --macro 10 --reference reference_0.pdb --grpname Protein --lag 2 --num_sims 4
#-----------------------------------------------------------------------------------------
