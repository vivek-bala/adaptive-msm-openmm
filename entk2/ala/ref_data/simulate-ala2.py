from __future__ import print_function
import sys
from simtk.openmm import app
import simtk.openmm as mm
from simtk import unit
from sys import stdout

class Tee(object):
    def __init__(self, name, mode):
        self.file = open(name, mode)
        self.stdout = sys.stdout
    def write(self, data):
        self.file.write(data)
        self.stdout.write(data)

pdb = app.PDBFile('ala2.pdb')
forcefield = app.ForceField('amber99sbildn.xml', 'amber99_obc.xml')

system = forcefield.createSystem(pdb.topology, nonbondedMethod=app.CutoffNonPeriodic,
    nonbondedCutoff=1.0*unit.nanometers, constraints=app.HBonds, rigidWater=True)
integrator = mm.LangevinIntegrator(300*unit.kelvin, 91.0/unit.picoseconds, 
    2.0*unit.femtoseconds)
integrator.setConstraintTolerance(0.00001)

platform = mm.Platform.getPlatformByName('CUDA')
properties = {'CudaPrecision': 'mixed'}
simulation = app.Simulation(pdb.topology, system, integrator, platform, properties)
simulation.context.setPositions(pdb.positions)

simulation.context.setVelocitiesToTemperature(300*unit.kelvin)
print('Equilibrating...')
simulation.step(2000)

initial_positions = simulation.context.getState(getPositions=True).getPositions()
interval = int(1*unit.picosecond / (2*unit.femtosecond))
n_steps = int(10*unit.nanosecond / (2*unit.femtosecond))

for i in range(10):
    print('Starting trajectory %d' % i)
    simulation.context.setPositions(initial_positions)
    simulation.reporters.append(app.DCDReporter('trajectory-%d.dcd' % i, interval))
    simulation.reporters.append(app.StateDataReporter(
        Tee('trajectory-%d.log' % i, 'w'), interval, step=True, time=True,
        potentialEnergy=True, temperature=True, progress=True,
        remainingTime=True, speed=True, totalSteps=n_steps, separator='\t'))
    
    simulation.step(n_steps)
    del simulation.reporters[:]

