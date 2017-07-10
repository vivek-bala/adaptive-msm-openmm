from __future__ import print_function
import sys
from simtk.openmm import app
import simtk.openmm as mm
from simtk import unit
from sys import stdout
import argparse
import mdtraj as md

class Tee(object):
    def __init__(self, name, mode):
        self.file = open(name, mode)
        self.stdout = sys.stdout
    def write(self, data):
        self.file.write(data)
        self.stdout.write(data)


if __name__ == '__main__':


    parser = argparse.ArgumentParser()
    parser.add_argument('--ns')

    args = parser.parse_args()

    pdb = md.load('ala2.pdb')
    forcefield = app.ForceField('amber99sbildn.xml', 'amber99_obc.xml')

    system = forcefield.createSystem(pdb.topology.to_openmm(), nonbondedMethod=app.CutoffNonPeriodic,
        nonbondedCutoff=1.0*unit.nanometers, constraints=app.HBonds, rigidWater=True)
    integrator = mm.LangevinIntegrator(300*unit.kelvin, 91.0/unit.picoseconds, 
        2.0*unit.femtoseconds)
    integrator.setConstraintTolerance(0.00001)

    platform = mm.Platform.getPlatformByName('Reference')
    #properties = {'CudaPrecision': 'mixed'}
    simulation = app.Simulation(pdb.topology.to_openmm(), system, integrator, platform)
    simulation.context.setPositions(pdb.xyz[0])

    simulation.context.setVelocitiesToTemperature(300*unit.kelvin)
    print('Equilibrating...')
    simulation.step(2000)

    #initial_positions = simulation.context.getState(getPositions=True).getPositions()
    interval = int(1*unit.picosecond / (2*unit.femtosecond))
    #interval = int(args.interval)
    n_steps = int(int(args.ns)*unit.nanosecond / (2*unit.femtosecond))
    #n_steps = int(args.n_steps)
    print (interval, n_steps)

    print('Starting trajectory')

    #simulation.context.setPositions(initial_positions)
    simulation.reporters.append(app.DCDReporter('trajectory.dcd', interval))
    simulation.reporters.append(app.StateDataReporter(
        Tee('trajectory.log', 'w'), interval, step=True, time=True,
        potentialEnergy=True, temperature=True, progress=True,
        remainingTime=True, speed=True, totalSteps=n_steps, separator='\t'))
    
    simulation.step(n_steps)
    del simulation.reporters[:]
