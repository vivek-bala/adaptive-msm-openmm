__author__    = "Vivek Balasubramanian <vivek.balasubramanian@rutgers.edu>"
__copyright__ = "Copyright 2016, http://radical.rutgers.edu"
__license__   = "MIT"

from radical.entk import PoE, EoP, AppManager, Kernel, ResourceHandle

from openmm import openmm_kernel
from echo import echo_kernel
from msm import msm_kernel

import argparse
import os
import re
import pprint
import glob


## USER PARS
ENSEMBLE_SIZE=25
PIPELINE_SIZE=2
INTERVAL=500  ## Picoseconds
N_STEPS=500000  
TOTAL_DUR = N_STEPS* INTERVAL / 1000   ## Nanoseconds

trial=2

class Test(PoE):

    def __init__(self, ensemble_size, pipeline_size, name):
        super(Test,self).__init__(ensemble_size, pipeline_size, name=name)

    def stage_1(self, instance):

        global ENSEMBLE_SIZE
        global INTERVAL
        global N_STEPS
            
        k1 = Kernel(name='openmm')
        k1.arguments = [ '--op_ind=%s'%instance, '--start_ind=1']
        k1.link_input_data = ['$SHARED/100-fs-peptide-400K.pdb','$SHARED/simulate-fs.py']
        return k1

    def stage_2(self, instance):

        k1 = Kernel(name="msm")
        k1.arguments = ['--lag=2', 
                        '--stride=10',
                        '--clusters=100',
                        '--components=4',
                        '--pdb=fs-peptide.pdb']

        k1.link_input_data = ['$SHARED/fs-peptide.pdb','$SHARED/analyze.py']

        for i in range(ENSEMBLE_SIZE):
            k1.link_input_data += ['$ITER_1_STAGE_1_TASK_%s/trajectory-%s.dcd > trajectory-%s.dcd'%(i+1, i+1, i+1)]

        k1.cores = 1

        return k1

if __name__ == '__main__':

    # Create pattern object with desired ensemble size, pipeline size
    pipe = Test(ensemble_size=[ENSEMBLE_SIZE,1], pipeline_size=PIPELINE_SIZE, name='1')

    # Create an application manager
    app = AppManager(name='MSM')

    # Register kernels to be used
    app.register_kernels(openmm_kernel)
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

        'xsede.stampede': {'cores': '64', 'username': 'vivek91','project': 'TG-MCB090174', 'queue': 'development', 'schema': 'gsissh'},
        'xsede.supermic': {'cores': '500', 'username': 'vivek91','project': 'TG-MCB090174', 'queue': 'workq', 'schema': 'gsissh'},
    }

    path=os.path.dirname(os.path.dirname(os.path.realpath(__file__)))

    # Create a resource handle for target machine
    res = ResourceHandle(resource=resource,
                cores=res_dict[resource]['cores'],
                username=res_dict[resource]['username'],
                project = res_dict[resource]['project'],
                queue= res_dict[resource]['queue'],
                walltime=60,
                database_url='mongodb://rp:rp@ds137749.mlab.com:37749/db_msm3',
                access_schema=res_dict[resource]['schema']
                )



    # Shared data
    res.shared_data = [ './100-fs-peptide-400K.pdb','./simulate-fs.py','./analyze.py','./fs-peptide.pdb']

    # Submit request for resources + wait till job becomes Active
    res.allocate(wait=True)

    try:

        # Run the first workload
        res.run(app)
        '''
        while(trial<=10):

            # Run the second workload
            ENSEMBLE_SIZE=100
            pipe2 = Test(ensemble_size=ENSEMBLE_SIZE+1, pipeline_size=PIPELINE_SIZE+1, name='{0}'.format(trial))

            app.add_workload(pipe2)

            USABLE_SIM_DATA = dict()
            USABLE_SIM_LIST = list()
            USABLE_SIM_ITER = [1 for x in range(1, ENSEMBLE_SIZE+1)]
            #USABLE_ANA_DATA = dict()
            MSM_FLAG=False
            MSM_DONE_FLAG=False

            ITER=[1 for x in range(1, ENSEMBLE_SIZE+2)]    

            res.run(app)
            trial+=1
        '''
    except Exception,ex:

        print 'Failed with error: ',ex

    finally:
        # Deallocate the resource
        res.deallocate()
