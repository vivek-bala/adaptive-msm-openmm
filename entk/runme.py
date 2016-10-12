__author__    = "Vivek Balasubramanian <vivek.balasubramanian@rutgers.edu>"
__copyright__ = "Copyright 2016, http://radical.rutgers.edu"
__license__   = "MIT"

from radical.entk import EoP, AppManager, Kernel, ResourceHandle

from grompp import grompp_kernel
from mdrun import mdrun_kernel
from traj_collect import traj_collect_kernel
from echo import echo_kernel
from msm_analysis import msm_kernel

import argparse
import os
import re


## USER PARS
ENSEMBLE_SIZE=4
PIPELINE_SIZE=3


CLUSTER_GEN=0
TERMINATE=False
TOTAL_TRAJ=0.0
RECLUSTER=1.0
RECLUSTER_NOW=True

ITER=[1 for x in range(1, ENSEMBLE_SIZE+2)]

class Test(EoP):

    def __init__(self, ensemble_size, pipeline_size):
        super(Test,self).__init__(ensemble_size, pipeline_size)

    def stage_1(self, instance):

        if instance <= ENSEMBLE_SIZE:
            k1 = Kernel(name='grompp')
            k1.arguments = [  
                                "--mdp=grompp.mdp",
                                "--conf=equil.gro",
                                "--top=topol.top",
                                "--out=topol.tpr"
                            ]        
            k1.cores=1

            k1.link_input_data = [
                                    '$SHARED/grompp.mdp',
                                    '$SHARED/equil{0}.gro > equil.gro'.format(instance-1),
                                    '$SHARED/topol.top'
                                ]

            return k1

        else:

            k1 = Kernel(name="echo")
            k1.arguments = ["--file=output.txt","--text=dummy"]
            k1.cores = 1

            return k1


    def stage_2(self, instance):

        if instance <= ENSEMBLE_SIZE:

            k2 = Kernel(name="mdrun")
            k2.arguments = [
                                "--tpr=topol.tpr",
                                "--rcon=0.7",
                                "--xtc=traj.xtc"
                            ]

            k2.cores=1
            k2.link_input_data = ['$STAGE_1_TASK_{0}/topol.tpr'.format(instance)]

            return k2

        else:

            k1 = Kernel(name="echo")
            k1.arguments = ["--file=output.txt","--text=dummy"]
            k1.cores = 1

            return k1


    def stage_3(self, instance):

        if instance <= ENSEMBLE_SIZE:

            k3 = Kernel(name="traj_collect")

            k3.arguments = [
                                '--xtc=traj.xtc',
                                '--system=Protein',
                                '--xtc_nopbc=traj.nopbc.xtc',
                                '--reference=reference.pdb',
                                '--lh5=file.lh5',
                                '--tpr=topol.tpr'
                ]


            k3.link_input_data = [
                                    '$STAGE_1_TASK_{0}/topol.tpr'.format(instance),
                                    '$STAGE_2_TASK_{0}/traj.xtc'.format(instance),
                                    '$SHARED/checktrajectory.py',
                                    '$SHARED/reference.pdb',
                                    '$SHARED/convert2lh5.py',
                                    '$SHARED/pre_analysis.py'
                                ]

            return k3

        else:

            k1 = Kernel(name="echo")
            k1.arguments = ["--file=output.txt","--text=dummy"]
            k1.cores = 1

            return k1


    def stage_4(self, instance):

        if instance <= ENSEMBLE_SIZE:

            k1 = Kernel(name="echo")
            k1.arguments = ["--file=output.txt","--text=dummy"]
            k1.cores = 1

            return k1

        else:

            m1 = Kernel(name="msm")
            m1.arguments = [
                                '--macro=10',
                                '--micro=100',
                                '--reference=reference_0.pdb',
                                '--grpname=Protein',
                                '--lag=2',
                                '--num_sims=20'
#                                '--ensembles=4'

                            ]
            m1.cores = 1

            m1.link_input_data = ['$SHARED/reference.pdb > reference_0.pdb',
                                '$SHARED/MSMproject.py']


            for i in range(1, ITER[instance-1]+1):
                for inst in range(1, ENSEMBLE_SIZE+1):

                    m1.link_input_data += [
                                        '$ITER_{1}_STAGE_3_TASK_{0}/traj.xtc > traj_{0}.xtc'.format(inst,i),
                                        '$ITER_{1}_STAGE_3_TASK_{0}/traj.nopbc.xtc > traj_{0}.nopbc.xtc'.format(inst,i),
                                        '$ITER_{1}_STAGE_3_TASK_{0}/topol.tpr > topol_{0}.tpr'.format(inst,i),
                                        '$ITER_{1}_STAGE_3_TASK_{0}/file.lh5 > file_{0}.lh5'.format(inst,i),
                                        '$ITER_{1}_STAGE_3_TASK_{0}/traj_info.txt > traj_info_{0}.txt'.format(inst,i)

                                    ]


            return m1

    def branch_4(self, instance):

        global TERMINATE
        global TOTAL_TRAJ
        global RECLUSTER
        global CLUSTER_GEN
        global RECLUSTER_NOW

        if instance <= ENSEMBLE_SIZE:

            if TERMINATE==False:
                self.set_next_stage(stage=1)  
            else:
                pass

        else:

            k = str(self.get_output(stage=4, instance=ENSEMBLE_SIZE+1))

            try:
                step = re.compile('^Total_ns=([0-9.]*)', re.MULTILINE)
                match = step.search(k)
                TOTAL_TRAJ+=float(match.group(1))
                print 'Traj: ',TOTAL_TRAJ, ' gen: ', CLUSTER_GEN
                print 'Diff:', (TOTAL_TRAJ - RECLUSTER*CLUSTER_GEN)
                if ((TOTAL_TRAJ - RECLUSTER*CLUSTER_GEN > RECLUSTER)or(RECLUSTER_NOW)):
                    print 'Total so far: ',TOTAL_TRAJ
                    ## Setup new simulations with new configurations
                    
                    ## Reiterate MSM in a while
                    self.set_next_stage(stage=4)
                    RECLUSTER_NOW=False
                else:

                    # Yaayyy !
                    print 'Reached convergence, terminate !'
                    TERMINATE=True

                CLUSTER_GEN+=1
            except:
                print 'Check-> Not enough simulations. Generate more !'
                self.set_next_stage(stage=4)




if __name__ == '__main__':

    # Create pattern object with desired ensemble size, pipeline size
    pipe = Test(ensemble_size=ENSEMBLE_SIZE+1, pipeline_size=PIPELINE_SIZE+1)

    # Create an application manager
    app = AppManager(name='MSM')

    # Register kernels to be used
    app.register_kernels(grompp_kernel)
    app.register_kernels(mdrun_kernel)
    app.register_kernels(traj_collect_kernel)
    app.register_kernels(echo_kernel)
    app.register_kernels(msm_kernel)


    # Add workload to the application manager
    app.add_workload(pipe)

    parser = argparse.ArgumentParser()
    parser.add_argument('--resource', help='target resource label')
    args = parser.parse_args()
    
    if args.resource != None:
        resource = args.resource
    else:
        resource = 'local.localhost'


    res_dict = {

        'xsede.stampede': {'cores': '16', 'username': 'vivek91','project': 'TG-MCB090174', 'queue': 'development', 'schema': 'gsissh'},
        'xsede.comet': {'cores': '24', 'username': 'vivek91', 'project': 'unc101', 'queue': 'compute', 'schema': 'gsissh'},
        'local.localhost': {'cores': '4', 'username': None, 'project': None, 'queue': None, 'schema': None},        

    }

    path=os.path.dirname(os.path.dirname(os.path.realpath(__file__)))

    # Create a resource handle for target machine
    res = ResourceHandle(resource=resource,
                cores=res_dict[resource]['cores'],
                username=res_dict[resource]['username'],
                project = res_dict[resource]['project'],
                queue= res_dict[resource]['queue'],
                walltime=20,
                database_url='mongodb://rp:rp@ds015335.mlab.com:15335/rp',
                access_schema=res_dict[resource]['schema']
                )



    # Shared data
    res.shared_data = [ '{0}/equil0.gro'.format(path),
                        '{0}/equil1.gro'.format(path),
                        '{0}/equil2.gro'.format(path),
                        '{0}/equil3.gro'.format(path),
                        './grompp.mdp',
                        '{0}/topol.top'.format(path),
                        '{0}/checktrajectory.py'.format(path),
                        './reference.pdb',
                        './convert2lh5.py',
                        './pre_analysis.py',
                        './MSMproject.py'
                    ]

    # Submit request for resources + wait till job becomes Active
    res.allocate(wait=True)

    # Run the given workload
    res.run(app)

    # Deallocate the resource
    res.deallocate()
