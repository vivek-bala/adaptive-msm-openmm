#!/bin/sh

## Setup GROMACS
export PATH=$PATH:/home/vivek/Research/modules/gromacs-5.1.3/build/bin

## Source location, cur location
src=/home/vivek/Research/repos/adaptive-msm
cur=$PWD

## Input parameters
system=Protein

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
## Tune and run simulations
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
    cd ..
    inst=$((inst+1))
done
#-----------------------------------------------------------------------------------------
