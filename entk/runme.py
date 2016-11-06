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
import pprint


## USER PARS
ENSEMBLE_SIZE=4
PIPELINE_SIZE=3
LAG_TIME=10
CLUSTER_GEN=1
TOTAL_TRAJ=0.0
RECLUSTER=1.0
#RECLUSTER_NOW=True


## INTERNAL PARS
USABLE_SIM_DATA = dict()
USABLE_SIM_LIST = list()
USABLE_SIM_ITER = [1 for x in range(1, ENSEMBLE_SIZE+1)]
USABLE_ANA_DATA = dict()
MSM_FLAG=False
MSM_DONE_FLAG=False

ITER=[1 for x in range(1, ENSEMBLE_SIZE+2)]

class Test(EoP):

    def __init__(self, ensemble_size, pipeline_size, name):
        super(Test,self).__init__(ensemble_size, pipeline_size, name=name)

    def stage_1(self, instance):

        global ENSEMBLE_SIZE
        global ITER
        global USABLE_ANA_DATA

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
                                    '$SHARED/topol.top'
                                ]

            if self._name == 'init':
                k1.link_input_data += ['$SHARED/equil{0}.gro > equil.gro'.format(instance-1)]
            else:
                k1.link_input_data += ['$PAT_init_ITER_{0}_STAGE_4_TASK_{1}/new_run_{2}.gro > equil.gro'.format(USABLE_ANA_DATA['iteration'],
                                                                                                                USABLE_ANA_DATA['instance'],
                                                                                                                instance-1
                                                                                                                )]
            
            return k1

        else:

            k1 = Kernel(name="echo")
            k1.arguments = ["--file=output.txt","--text=st1_msm"]
            k1.cores = 1
            return k1


    def stage_2(self, instance):

        global ENSEMBLE_SIZE
        global ITER

        if instance <= ENSEMBLE_SIZE:

            k2 = Kernel(name="mdrun")
            k2.arguments = [
                                "--tpr=topol.tpr",
                                "--rcon=0.7",
                                "--xtc=traj.xtc"
                            ]

            k2.cores=1
            k2.link_input_data = ['$ITER_{1}_STAGE_1_TASK_{0}/topol.tpr'.format(instance, ITER[instance-1])]
            print 'HERE s2, iter:{0}, inst:{1} '.format(ITER[instance-1], instance)

            #pprint.pprint(self._pattern_dict)

            #print str(self.get_output(stage=1, instance=instance))

            return k2

        else:

            k1 = Kernel(name="echo")
            k1.arguments = ["--file=output.txt","--text=st2_msm"]
            k1.cores = 1
            return k1


    def stage_3(self, instance):

        global ENSEMBLE_SIZE
        global ITER

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
                                    '$ITER_{1}_STAGE_1_TASK_{0}/topol.tpr'.format(instance, ITER[instance-1]),
                                    '$ITER_{1}_STAGE_2_TASK_{0}/traj.xtc'.format(instance, ITER[instance-1]),
                                    '$SHARED/checktrajectory.py',
                                    '$SHARED/reference.pdb',
                                    '$SHARED/convert2lh5.py',
                                    '$SHARED/pre_analysis.py'
                                ]

            return k3

        else:

            k1 = Kernel(name="echo")
            k1.arguments = ["--file=output.txt","--text=st3_msm"]
            k1.cores = 1
            return k1

    def branch_3(self, instance):


        # Logging vars
        global USABLE_SIM_DATA
        global USABLE_SIM_LIST
        global ITER
        global MSM_FLAG

        # Decision making vars
        global TOTAL_TRAJ
        global RECLUSTER
        global CLUSTER_GEN
        global RECLUSTER_NOW

        if ((instance <= ENSEMBLE_SIZE)and(MSM_FLAG==False)and(MSM_DONE_FLAG==False)):


            try:

                # Try-catch unreliable hack to get thru race condition
                k = str(self.get_output(stage=3, instance=instance))
                step = re.compile('^Total_ns=([0-9.]*)', re.MULTILINE)
                match = step.search(k)
                TOTAL_TRAJ+=float(match.group(1))

                if float(match.group(1)) > 0.0:

                    # Note the required sims
                    USABLE_SIM_DATA = {'instance': instance, 'iter': ITER[instance-1]}
                    USABLE_SIM_LIST.append(USABLE_SIM_DATA)
            except:
                pass

            print 'Total traj: ', TOTAL_TRAJ
            print 'Diff: ', (TOTAL_TRAJ - RECLUSTER*CLUSTER_GEN),', RECLUSTER: ', RECLUSTER


            #if ((TOTAL_TRAJ - RECLUSTER*CLUSTER_GEN > RECLUSTER)or(RECLUSTER_NOW)):
            if ((TOTAL_TRAJ - RECLUSTER*CLUSTER_GEN > RECLUSTER)):

                print 'Total traj: ', TOTAL_TRAJ
                print 'Diff: ', (TOTAL_TRAJ - RECLUSTER*CLUSTER_GEN),', RECLUSTER: ', RECLUSTER
                print 'Instance: {0}, Iter: {1}, Dict:'.format(instance,ITER[instance-1],self._pattern_dict)

                # Cancel currently runnings sim tasks
                #if RECLUSTER_NOW != True:
                #self.cancel_all_tasks = True

                # Run MSM
                print 'Going to MSM'
                self.set_next_stage(stage=4)
                MSM_FLAG=True

                # Reassignments
                CLUSTER_GEN+=1
                #RECLUSTER_NOW=False

            else:

                print 'Skip.Instance: {0}, Iter: {1}'.format(instance,ITER[instance-1],self._pattern_dict)

                self.set_next_stage(stage=1)            
                ITER[instance-1] += 1

        else:

            pass


    def stage_4(self, instance):

        global MSM_FLAG
        global ENSEMBLE_SIZE
        global USABLE_SIM_LIST
        global USABLE_ANA_DATA
        global MSM_DONE_FLAG

        if instance <= ENSEMBLE_SIZE:

            print 'Instance {0} should end'.format(instance)
            k1 = Kernel(name="echo")
            k1.arguments = ["--file=output.txt","--text=st4_inst"]
            k1.cores = 1
            return k1


        if MSM_FLAG == True:
            print 'This is MSM'

            USABLE_ANA_DATA['instance'] = instance 
            USABLE_ANA_DATA['iteration'] = ITER[instance-1]

            m1 = Kernel(name="msm")
            m1.arguments = [
                                '--macro=10',
                                '--micro=100',
                                '--reference=reference_0.pdb',
                                '--grpname=Protein',
                                '--lag=2',
                                '--num_sims=20'
                            ]
            m1.cores = 1

            m1.link_input_data = ['$SHARED/reference.pdb > reference_0.pdb',
                                '$SHARED/MSMproject.py']

            for item in USABLE_SIM_LIST:

                inst = item['instance']
                i = item['iter']

                index = str(inst)+str(i)

                m1.link_input_data += [
                                        '$ITER_{0}_STAGE_3_TASK_{1}/traj.xtc > traj_{2}.xtc'.format(i, inst, index),
                                        '$ITER_{0}_STAGE_3_TASK_{1}/traj.nopbc.xtc > traj_{2}.nopbc.xtc'.format(i, inst, index),
                                        '$ITER_{0}_STAGE_3_TASK_{1}/topol.tpr > topol_{2}.tpr'.format(i, inst, index),
                                        '$ITER_{0}_STAGE_3_TASK_{1}/file.lh5 > file_{2}.lh5'.format(i, inst, index),
                                        '$ITER_{0}_STAGE_3_TASK_{1}/traj_info.txt > traj_info_{2}.txt'.format(i, inst, index)
                                        ]

                m1.copy_output_data = []
                for i in range(0,200):
                    m1.copy_output_data += ['new_run_{0}.gro > $SHARED/new_run_{0}.gro'.format(i)]


            MSM_FLAG=False
            MSM_DONE_FLAG=True

            return m1

        else:

            print 'This is dummy MSM'

            k1 = Kernel(name="echo")
            k1.arguments = ["--file=output.txt","--text=st4_msm"]
            k1.cores = 1
            return k1


    def branch_4(self, instance):

        global ITER
        global ENSEMBLE_SIZE
        global MSM_FLAG
        global MSM_DONE_FLAG


        try:
            k = str(self.get_output(stage=4, instance=instance))
            step = re.compile('^Total_ns=([0-9.]*)', re.MULTILINE)
            match = step.search(k)

            #if float(match.group(1)) > 0.0:
            #    ENSEMBLE_SIZE = 200
            #    self._ensemble_size = ENSEMBLE_SIZE+1

        except:
            pass
            
        if ((instance>ENSEMBLE_SIZE)and(MSM_DONE_FLAG == False)):
            self.set_next_stage(stage=1)
            ITER[instance-1] += 1



if __name__ == '__main__':

    # Create pattern object with desired ensemble size, pipeline size
    pipe = Test(ensemble_size=ENSEMBLE_SIZE+1, pipeline_size=PIPELINE_SIZE+1, name='init')

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

    # Run the first workload
    res.run(app)


    # Run the second workload
    #pipe2 = Test(ensemble_size=200+1, pipeline_size=PIPELINE_SIZE+1, name='sec')

    #app.add_workload(pipe2)

    #res.run(app)

    # Deallocate the resource
    res.deallocate()
