from __future__ import print_function
import sys
import mdtraj as md
from simtk.openmm import app
import simtk.openmm as mm
from simtk import unit
import argparse

class Tee(object):
    def __init__(self, name, mode):
        self.file = open(name, mode)
        self.stdout = sys.stdout
    def write(self, data):
        self.file.write(data)
        self.stdout.write(data)


if __name__ == '__main__':

    if len(sys.argv) != 3:
        print('usage %s <trajectory index (for output file)> <model index of starting conformation>')
        exit(1)

    pdb = md.load('100-fs-peptide-400K.pdb')
    forcefield = app.ForceField('amber99sbildn.xml', 'amber99_obc.xml')

    system = forcefield.createSystem(pdb.topology.to_openmm(), nonbondedMethod=app.CutoffNonPeriodic,
        nonbondedCutoff=1.0*unit.nanometers, constraints=app.HBonds)
    integrator = mm.LangevinIntegrator(300*unit.kelvin, 91.0/unit.picoseconds, 
        2.0*unit.femtoseconds)
    integrator.setConstraintTolerance(0.00001)

    platform = mm.Platform.getPlatformByName('CPU')
    #properties = {'CudaPrecision': 'mixed', 'CudaDeviceIndex': sys.argv[1]}
    simulation = app.Simulation(pdb.topology.to_openmm(), system, integrator, platform)
    simulation.context.setPositions(pdb.xyz[int(sys.argv[2])])
    
    simulation.context.setVelocitiesToTemperature(300*unit.kelvin)
    
    nsteps = int((1*unit.nanoseconds) / (2*unit.femtoseconds))
    interval = int((10*unit.picoseconds) / (2*unit.femtoseconds))
    
    simulation.reporters.append(app.StateDataReporter(open('trajectory-%s.log' % sys.argv[1], 'w', 0),
        interval, step=True, time=True, progress=True,
        potentialEnergy=True, temperature=True, remainingTime=True,
        speed=True, totalSteps=nsteps, separator='\t'))
    
    # equilibrate
    simulation.step(int(100*unit.picoseconds / (2*unit.femtoseconds)))
    
    # now add the trajectory reporter.
    simulation.reporters.append(app.DCDReporter('trajectory-%s.dcd' % sys.argv[1], interval))
    simulation.step(nsteps)
