# This file is part of Copernicus
# http://www.copernicus-computing.org/
# 
# Copyright (C) 2011, Sander Pronk, Iman Pouya, Erik Lindahl, and others.
#
# This program is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License version 2 as published 
# by the Free Software Foundation
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License along
# with this program; if not, write to the Free Software Foundation, Inc.,
# 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.



import os
import sys
import subprocess
import re
import logging

import traceback
try:
    from cStringIO import StringIO
except ImportError:
    from StringIO import StringIO


from scipy.stats import halfnorm
import scipy 
import scipy.sparse
import random
from numpy import where
from numpy import array,argmax
import numpy
#import random

#import matplotlib
#Use a non GUI backend for matplotlib
#matplotlib.use('Agg')
#import matplotlib.pyplot as plt


# Make sure msmbuilder is in the PYTHONPATH
from msmbuilder.CopernicusProject import *
import msmbuilder.MSMLib
import msmbuilder.Serializer
import msmbuilder.Trajectory
#import DataFile

import argparse
import glob

log=logging.getLogger(__name__)


class TrajData(object):
    """Information about a trajectory"""
    def __init__(self, lh5, xtc, xtc_nopbc, tpr, dt, frames):
        self.lh5=lh5
        self.xtc=xtc
        self.xtc_nopbc=xtc_nopbc
        self.tpr=tpr
        self.dt=dt
        self.frames=frames


class MSMProject(object):
    ''' MSM specific project data '''

    # the name of the top-level element for this project data type
    elementName = ""

    def __init__(self, microstates, macrostates, reference, grpname, lag_time, num_sims, index_file=None):
        ''' Initialize MSMProject '''
        
        self.num_micro     = int(microstates)
        if self.num_micro <= 1:
            sys.stderr("Error: num_micro=%d: invalid number\n"%self.num_micro)
        self.num_macro     = int(macrostates)
        if self.num_macro <= 1:
            sys.stderr("Error: num_macro=%d: invalid number\n"%self.num_macro)
        self.ref_conf      = reference
        self.grpname       = grpname
        #self.num_sim       = inp.getInput('num_to_start')
        self.lag_time      = int(lag_time)

        if self.lag_time is not None and self.lag_time <= 0:
            sys.stderr("Error: lag_time=%g: invalid number\n"%self.lag_time)

        # Set index file parameter
        self.ndx = index_file

        # The msm-builder project
        self.Proj = None

        # The assignments from msm-builder
        self.assignments = None

        # The transition count matrix from msm-builder
        self.T = None

        # The desired lag time
        self.max_time = None

        # Sims to start per round
        self.num_to_start = int(num_sims)

        # handle trajectories
        self.avgtime=0.
        self.filelist=[]
        self.trajData=dict()

        self._flag=True
        #ta=self.inp.getInput('trajectories')
        i=1

        for file in glob.glob('*.lh5'):
            index = file.split('.')[0].strip().split('_')[1].strip()
            lh5='file_{0}.lh5'.format(index)
            xtc='traj_{0}.xtc'.format(index)
            xtc_nopbc='traj_{0}.nopbc.xtc'.format(index)
            tpr='topol_{0}.tpr'.format(index)
            self.tpr = tpr

            try:
                with open('traj_info_{0}.txt'.format(index),'r') as f:
                    lines = f.readlines()
                    ns = float(lines[0].strip().split('=')[1].strip())
                    dt = float(lines[1].strip().split('=')[1].strip())
                    frames = int(lines[2].strip().split('=')[1].strip())
            except:
                self._flag=False

            self.filelist.append([lh5])
            self.trajData[lh5]=TrajData(lh5, xtc, xtc_nopbc, tpr, dt, frames)
            self.avgtime += dt * (frames-1)/1000.
            i+=1
        self.avgtime /= i
        sys.stderr.write("Average trajectory time: %g ns\n"%(self.avgtime))
        sys.stderr.write("filelist size=%d.\n"%(len(self.filelist)))

        random.seed()  

        if len(self.filelist)==0:
            self._flag=False


    @property
    def flag(self):
        return self._flag

    def getNewSimTime(self):
        ''' Compute a new simulation time from a half-normal distribution '''
     
        # Extend to 400 ns (hardcoded for villin)
        new_length = 400000
        
        r = random.random()

        if(r>0.9):
            nst = int(new_length/self.dt)
        else:
            nst = 25000000
    
        return nst


    def createMicroStates(self):
        ''' Build a micro-state MSM '''

        # Create the msm project from the reference conformation
        #TODO IMAN provide weighting here
        Proj = CreateCopernicusProject(self.ref_conf, self.filelist)
        self.Proj = Proj
        C1   = Conformation.Conformation.LoadFromPDB(self.ref_conf)
        
        # Automate the clustering to only CA or backbone atoms
        # TODO: fix this
        a = C1["AtomNames"]
        AtomIndices=where((a=="N") | (a=="C") | (a=="CA") | (a=="O"))[0]
        
        sys.stderr.write("Cluster project.\n")
        # Do msm-stuff
        GenF = os.path.join('Data','Gens.nopbc.h5')
        AssF = os.path.join('Data','Ass.nopbc.h5')
        AssFTrimmed = os.path.join('Data','Assignment-trimmed.nopbc.h5')
        RmsF = os.path.join('Data','RMSD.nopbc.h5')
        
        Generators = Proj.ClusterProject(AtomIndices=AtomIndices,
                                         NumGen=self.num_micro,Stride=30)
        sys.stderr.write("Assign project.\n")
        Assignments,RMSD,WhichTrajs = Proj.AssignProject(Generators,
                                                       AtomIndices=AtomIndices)
        if os.path.exists(GenF):
            os.remove(GenF)
        Generators.SaveToHDF(GenF)
        if os.path.exists(AssF):
            os.remove(AssF)
        msmbuilder.Serializer.SaveData(AssF,Assignments)
        if os.path.exists(RmsF):
            os.remove(RmsF)
        msmbuilder.Serializer.SaveData(RmsF,RMSD)
        
        print "Trim data.\n"
        # Trim data
        Counts = msmbuilder.MSMLib.GetCountMatrixFromAssignments(Assignments,
                                                       self.num_micro,
                                                       LagTime=1,
                                                       Slide=True)
                
        # Get the most populated state
        print "Get the most populated state.\n"
        X0       = array((Counts+Counts.transpose()).sum(0)).flatten()
        X0       = X0/sum(X0)
        MaxState = argmax(X0)


        ## Calculate only times up to at maximum half the
        ## length of an individual trajectory
        max_time = self.avgtime/2.
        print 'avgtime: ',self.avgtime
        #max_time = ((self.dt * self.nstep / 1000)*0.5)
        ## SP this is almost certainly wrong:
        #if max_time > 1:
        #    max_time=int(max_time)
        #else:
        #    max_time=2
        ###max_time = 300 # hard-coded for villin

        self.max_time = max_time
        
        # More trimming
        # PK want ErgodicTrim instead of EnforceMetastability
        # This is from BuildMSM script
        print "More trimming...\n"
        CountsAfterTrimming,Mapping=msmbuilder.MSMLib.ErgodicTrim(Counts)
        msmbuilder.MSMLib.ApplyMappingToAssignments(Assignments,Mapping)
        ReversibleCounts = msmbuilder.MSMLib.IterativeDetailedBalance(CountsAfterTrimming,Prior=0)
        TC = msmbuilder.MSMLib.EstimateTransitionMatrix(ReversibleCounts)
        Populations=numpy.array(ReversibleCounts.sum(0)).flatten()
        Populations/=Populations.sum()

        self.assignments=Assignments
        self.T=TC

        print 'assign: ',Assignments.flatten()

        NumStates=max(Assignments.flatten())+1
        print "New number of states=%d\n"%NumStates
        if os.path.exists(AssFTrimmed):
            os.remove(AssFTrimmed)
        msmbuilder.Serializer.SaveData(AssFTrimmed,Assignments)
        
        print "Calculating implied time scales..\n"
        # Calculate the implied time-scales
        time = numpy.arange(1,max_time+1,1)


        TS = msmbuilder.MSMLib.GetImpliedTimescales(AssFTrimmed,NumStates,time,NumImpliedTimes=len(time)+1)
        print "TS=%s, time=%s\n"%(str(TS), time)
        '''
        try:
            plt.scatter(TS[:,0],TS[:,1])
            plt.title('Lag times versus implied time scale')
            plt.xlabel('Lag Time (assignment-steps)')
            plt.ylabel('Implied Timescale (ps)')
            plt.yscale('log')
            timescalefn= 'msm_timescales.png'
            sys.stderr.write('Writing timescale plot to %s'%timescalefn)
            try:
                plt.savefig(timescalefn)
            except:
                fo=StringIO()
                traceback.print_exception(sys.exc_info()[0], 
                                          sys.exc_info()[1],
                                          sys.exc_info()[2], file=fo)
                errmsg="Run error generating timescale plot: %s\n"%(fo.
                                                                    getvalue())
                sys.stderr.write(errmsg)
            plt.close()
            self.out.setOut('timescales', FileValue(timescalefn))
        except ValueError as e:
            fo=StringIO()
            traceback.print_exception(sys.exc_info()[0], sys.exc_info()[1],
                                      sys.exc_info()[2], file=fo)
            errmsg="Run error generating timescale plot: %s\n"%(fo.getvalue())
            sys.stderr.write(errmsg)

        '''
        # Get random confs from each state
        sys.stderr.write("Getting random configuration from each state..\n")
        RandomConfs = Proj.GetRandomConfsFromEachState(Assignments,NumStates,1,
                                                       JustGetIndices=True)
        
        # Compute the MaxState with the new assignments (ie. after trimming)
        sys.stderr.write("Computing MaxState.\n")
        Counts=msmbuilder.MSMLib.GetCountMatrixFromAssignments(Assignments,
                                                               NumStates,
                                                               LagTime=1,
                                                               Slide=True)
        X0=array((Counts+Counts.transpose()).sum(0)).flatten()
        X0=X0/sum(X0)
        MaxState=argmax(X0)

        # Create a tpr-file for trjconv with -pbc mol
        #sys.stderr.write("making randomconfs.\n")
        #try:
        #    os.mkdir('RandomConfs')
        #except:
        #    pass
        # we need a tpr file to be able to trjconv random confs later
        #proc = subprocess.Popen(["grompp","-f","%s"%self.mdpfile,
        #                          "-c","%s"%self.grofile[0],
        #                          "-p", "%s"%self.topfile,"-o",
        #                          "%s"%os.path.join(self.inp.getOutputDir(),
        #                                            'topol.tpr')],
        #                      stdin=None,stdout=sys.stdout, stderr=sys.stdout)
        #proc.communicate(None)

        # we pick one of the tpr files.
        self.tprfile=self.tpr
        

        # Set a flag to indicate if we have written the maxstate.pdb-file
        have_maxstate=0
        
        for i in xrange(NumStates):
            traj_num    = RandomConfs[i][0][0]
            frame_nr    = RandomConfs[i][0][1]
            lh5name     = Proj.GetTrajFilename(traj_num)            
            #sys.stderr.write("trajectory name=%s\n"%lh5name)
            trajdata    = self.trajData[lh5name]
            trajname    = trajdata.xtc
            #trajname    = trajname.replace('.nopbc.lh5','.xtc')
            time        = frame_nr * trajdata.dt #* self.nstxtcout

            #if(i<10*self.num_to_start):
                #proc = subprocess.Popen(["trjconv","-f","%s"%trajname,"-s","%s"%os.path.join(self.inp.getOutputDir(),'topol.tpr'),"-o",os.path.join(self.inp.getOutputDir(),'micro%d.gro'%i),"-pbc","mol","-dump","%d"%time], stdin=subprocess.PIPE, stdout=sys.stdout, stderr=sys.stderr)
                
                #proc.communicate("0")

            # Write out a pdb of the most populated state
            if(i==MaxState and have_maxstate==0):
                maxstatefn='maxstate_{0}.pdb'.format(i)
                sys.stderr.write("writing out pdb of most populated state.\n")
                args1 = 'gmx trjconv'
                args2 = ["-f", trajname, "-s", self.tprfile,
                         "-o", maxstatefn, "-pbc", "mol", "-dump", "%d" % time]

                for item in args2:
                    args1 = args1 + ' ' + item

                if self.ndx is not None:
                    args1 += "-n {0}".format(self.ndx)

                os.system('echo {1} | {0}'.format(args1, self.grpname))
                #proc = subprocess.Popen(args, stdin=subprocess.PIPE, stdout=sys.stdout, stderr=sys.stderr)
                #[std, err] = proc.communicate(self.grpname)

                #print std
                #print err

                #self.out.setOut('maxstate', FileValue(maxstatefn))
                have_maxstate=1

        # now evenly sample configurations and put them in the array
        # newRuns. If we're later assigning macrosates, we'll overwrite them
        # with adaptive sampling configurations
        self.newRuns=[]
        for j in xrange(self.num_to_start*self.num_macro):
            # pick a cluster at random:
            i=random.random()*int(NumStates)
            traj_num    = int(RandomConfs[i][0][0])
            frame_nr    = int(RandomConfs[i][0][1])
            lh5name     = Proj.GetTrajFilename(traj_num)            
            trajdata    = self.trajData[lh5name]
            trajname    = trajdata.xtc
            time        = frame_nr * trajdata.dt 
            #maxstatefn=os.path.join(self.inp.getOutputDir(), '.conf')
            outfn= 'new_run_%d.gro'%(j)
            args1 = 'gmx trjconv'
            args2 = ["-f", "%s"%trajname, "-s", self.tprfile, 
                     "-o", outfn, "-pbc", "mol", "-dump", "%d" % time]

            for item in args2:
                    args1 = args1 + ' ' + item

            sys.stderr.write("writing out new run %s .\n"%outfn)
            #proc = subprocess.Popen(args, stdin=subprocess.PIPE, stdout=sys.stdout, stderr=sys.stderr)
            os.system('echo 0 | {0}'.format(args1))
            #proc.communicate('0')
            self.newRuns.append(outfn)

        print 'done with micro-states'



        #os.remove('mdout.mdp')

        # Make a plot of the rmsd vs rel. population (rmsd)
#        NumConfsPerState=1
#        RandomConfs = Proj.GetRandomConfsFromEachState(Assignments,NumStates,NumConfsPerState,JustGetIndices=False)
#        Allatoms=RandomConfs["Atoms"]

#        CA=intersect1d(AtomRange,where(Allatoms=="CA")[0])
#        rmsd=RandomConfs.CalcRMSD(C1,CA,CA).reshape((NumStates,NumConfsPerState)).mean(1)
 
 #       NumEigen=NumStates/100
 #       EigAns=msmbuilder.MSMLib.GetEigenvectors(T,NumEigen);
 #       Populations=EigAns[1][:,0]

 #       plt.plot(rmsd,-log(Populations),'o')
 #       plt.title("Free Energy Versus RMSD [nm]")
 #       plt.ylabel("Free Energy")
 #       plt.xlabel("RMSD [nm]")
 #       plt.savefig(os.path.join('cpc-data','msm_fe.png'))
 #       plt.close()


    def createMacroStates(self):
        ''' Build a macro-state MSM '''
        # Again we redirect output
        #stdoutfn=os.path.join(self.inp.getOutputDir(), 'msm_stdout_macro.txt')
        #stderrfn=os.path.join(self.inp.getOutputDir(), 'msm_stderr_macro.txt')
        
        #old_stdout = sys.stdout
        #sys.stdout=open(stdoutfn,'w')
        #old_stderr = sys.stderr
        #sys.stderr=open(stderrfn,'w')  
        
        Map         = msmbuilder.MSMLib.PCCA(self.T,self.num_macro)
        Assignments = self.assignments
        Assignments = Map[Assignments]
        NumStates = max(Assignments.flatten())+1

        sys.stderr.write("Calculating macrostates with lag time %g.\n"%self.lag_time)

        # Now repeat any calculations with the new assignments
        Counts = msmbuilder.MSMLib.GetCountMatrixFromAssignments(Assignments, self.num_macro, LagTime=self.lag_time, Slide=True)
        
        #PK want reversible MLE estimator again here
        sys.stderr.write("Recalculating assignments & trimming again.\n")
        CountsAfterTrimming,Mapping=msmbuilder.MSMLib.ErgodicTrim(Counts)
        msmbuilder.MSMLib.ApplyMappingToAssignments(Assignments,Mapping)
        ReversibleCounts = msmbuilder.MSMLib.IterativeDetailedBalance(
                                                         CountsAfterTrimming,
                                                         Prior=0)
        TC = msmbuilder.MSMLib.EstimateTransitionMatrix(ReversibleCounts)
        Populations=numpy.array(ReversibleCounts.sum(0)).flatten()
        Populations/=Populations.sum()

        # Again, get the most populated state
        X0       = array((Counts+Counts.transpose()).sum(0)).flatten()
        X0       = X0/sum(X0)
        MaxState = argmax(X0)

        tcoutf= "tc.dat"
        if scipy.sparse.issparse(TC):
            scipy.savetxt(tcoutf, TC.todense())
        else:
            numpy.savetxt(tcoutf, TC, fmt="%12.6g" )
        #self.out.setOut('macro_transition_counts', FileValue(tcoutf))

        woutf= "weights.dat"
        numpy.savetxt(woutf, X0, fmt="%12.6g" )
        #self.out.setOut('macro_weights', FileValue(woutf))

        with open('sim_details.txt','w') as f:
            f.write("macro_transition_counts = {0}".format(tcoutf))
            f.write("macro_weights = {0}".format(woutf))

       
        # Do adaptive sampling on the macrostates
        nstates=int(self.num_macro*self.num_to_start)
        sys.stderr.write("Adaptive sampling to %d=%d*%d states.\n"%(nstates, self.num_macro, self.num_to_start))
        Proj = self.Proj

        StartStates = Proj.AdaptiveSampling(Counts.toarray(),nstates)

        #print StartStates

        #PK note JustGetIndices gives indices into original conformations
        RandomConfs = Proj.GetRandomConfsFromEachState(Assignments,NumStates,1, JustGetIndices=True)
       
        self.newRuns=[]
        self.macroConfs=[]
        for k,v in StartStates.items():
            num_started = 0
            for i in xrange(NumStates):
                if i==k:
                    trajnum  = int(RandomConfs[i][0][0])
                    frame_nr = int(RandomConfs[i][0][1])
                    lh5name  = Proj.GetTrajFilename(trajnum)            
                    trajdata    = self.trajData[lh5name]
                    trajname    = trajdata.xtc
                    time        = frame_nr * trajdata.dt #* self.nstxtcout
                    #time     = frame_nr * self.dt *self.nstxtcout
                    #trajname = Proj.GetTrajFilename(trajnum)
                    #trajname = trajname.replace('.nopbc.lh5','.xtc')

                    first=True
                    # Use trjconv to write new starting confs
                    while(num_started < self.num_to_start):
                        sys.stderr.write("Writing new start confs.\n")
                        outfn=os.path.join('macro%d-%d.gro'%(i,num_started))

                        args1 = 'gmx trjconv'
                        args2 = ["-f", "%s" % trajname, "-s", self.tprfile,"-o", outfn, "-pbc", "mol", "-dump", "%d" % time]

                        for item in args2:
                            args1 = args1 + ' ' + item


                        os.system('echo 0 | {0}'.format(args1))
                        num_started = num_started + 1
                        self.newRuns.append(outfn)
                        if first:
                            self.macroConfs.append(outfn)
                            first=False

        # now set the macro state outputs:
        i=0
        for fname in self.macroConfs:
            with open('sim_details.txt','a') as f:
                f.write("macro_conf_{1} = {0}".format(fname,i))
            i+=1

        print 'done with macro-states'



    def calc_total_ns(self):

        total_traj_ns = 0.0

        for filename in glob.glob('traj_info_*.txt'):

            with open(filename,'r') as f:
                lines = f.readlines()
                ns = float(lines[0].strip().split('=')[1].strip())
                dt = float(lines[1].strip().split('=')[1].strip())
                frames = int(lines[2].strip().split('=')[1].strip())

                ns=dt*(frames-1)/1000.
                    
                total_traj_ns += ns

        print 'Total_ns={0}'.format(total_traj_ns)
        
if __name__ == '__main__':

    parser = argparse.ArgumentParser()
    parser.add_argument("--micro", help="number of micro states")
    parser.add_argument("--macro", help="number of macro states")
    parser.add_argument("--reference", help="reference filename")
    parser.add_argument("--grpname", help="groupname")
    parser.add_argument("--lag", help="lag time")
    parser.add_argument("--num_sims", help="number of simulations per state")
    #parser.add_argument("--ensembles", help="number of min simulations")
    args = parser.parse_args()


    #if ((len(glob.glob('*.tpr'))>=int(args.ensembles)) and (len(glob.glob('*.nopbc.xtc'))>=int(args.ensembles)) 
    #    and (len(glob.glob('*.xtc'))>=int(args.ensembles)) and (len(glob.glob('*.lh5'))>=int(args.ensembles))):


    msmproject = MSMProject(microstates=args.micro, 
                            macrostates=args.macro, 
                            reference=args.reference, 
                            grpname=args.grpname, 
                            lag_time=args.lag, 
                            num_sims=args.num_sims)


    if msmproject.flag != False:
        # Build the microstates
        msmproject.createMicroStates()

        # Build the macrostates
        msmproject.createMacroStates()

        # Compute total traj ns
        msmproject.calc_total_ns()

    else:

        print 'Error ! Probably a read error of one of the files !'
