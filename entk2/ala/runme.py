__author__    = "Vivek Balasubramanian <vivek.balasubramanian@rutgers.edu>"
__copyright__ = "Copyright 2016, http://radical.rutgers.edu"
__license__   = "MIT"

from radical.entk import PoE, EoP, AppManager, Kernel, ResourceHandle

from openmm import openmm_kernel
from msm import msm_kernel

import argparse
import os
import re
import pprint
import glob


## USER PARS
ENSEMBLE_SIZE=int(os.environ.get('ENSEMBLE_SIZE',None))
PIPELINE_SIZE=2
NS=int(os.environ.get('NS',None))
ITER=1
TOTAL_ITERS=int(os.environ.get('TOTAL_ITERS',None))

class Test(PoE):

    def __init__(self, ensemble_size, pipeline_size, name, iterations):
        super(Test,self).__init__(ensemble_size, pipeline_size, name=name, iterations=iterations)

    def stage_1(self, instance):

        global ENSEMBLE_SIZE
        global INTERVAL
        global N_STEPS, ITER
            
        k1 = Kernel(name='openmm')
        k1.arguments = [ '--ns=%s'%NS]
        if ITER==1:
            k1.link_input_data = ['$SHARED/ala2.pdb','$SHARED/simulate.py']
        else:
            k1.link_input_data = ['$ITER_%s_STAGE_2_TASK_1/ala2-%s.pdb > ala2.pdb'%(ITER-1, instance-1), '$SHARED/simulate.py']
        return k1

    def stage_2(self, instance):


        global ITER, NS, ENSEMBLE_SIZE, TOTAL_ITERS
        k1 = Kernel(name="msm")
        k1.arguments = ['--lag=2', 
                        '--stride=10',
                        '--clusters=100',
                        '--components=4',
                        '--pdb=ala2.pdb']

        k1.link_input_data = ['$SHARED/ala2.pdb','$SHARED/analyze.py']

        for i in range(ITER):
            for j in range(ENSEMBLE_SIZE):
                k1.link_input_data += ['$ITER_%s_STAGE_1_TASK_%s/trajectory.dcd > trajectory-%s_%s.dcd'%(i+1,j+1,i,j)]

        k1.cores = 1

        k1.download_output_data = ['microstate_info.txt > dur-%s-ensemble-%s-iters-%s/microstate_info-%s.txt'%(NS, ENSEMBLE_SIZE, TOTAL_ITERS, ITER),
                                   'macrostate_info.txt > dur-%s-ensemble-%s-iters-%s/macrostate_info-%s.txt'%(NS, ENSEMBLE_SIZE, TOTAL_ITERS, ITER)]

        ITER+=1

        return k1
    

if __name__ == '__main__':

    # Create pattern object with desired ensemble size, pipeline size
    pipe = Test(ensemble_size=[ENSEMBLE_SIZE,1], pipeline_size=PIPELINE_SIZE, name='1', iterations=TOTAL_ITERS)

    # Create an application manager
    app = AppManager(name='MSM')

    # Register kernels to be used
    app.register_kernels(openmm_kernel)
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

        'xsede.stampede': {'cores': ENSEMBLE_SIZE, 'username': 'vivek91','project': 'TG-MCB090174', 'queue': 'development', 'schema': 'gsissh'},
        'xsede.supermic': {'cores': ENSEMBLE_SIZE, 'username': 'vivek91','project': 'TG-MCB090174', 'queue': 'workq', 'schema': 'gsissh'},
    }

    path=os.path.dirname(os.path.dirname(os.path.realpath(__file__)))

    # Create a resource handle for target machine
    res = ResourceHandle(resource=resource,
                cores=res_dict[resource]['cores'],
                username=res_dict[resource]['username'],
                project = res_dict[resource]['project'],
                queue= res_dict[resource]['queue'],
                walltime=480,        # Roughly 1.2 mins/ns for this system
                database_url='mongodb://rp:rp@ds137749.mlab.com:37749/db_msm3',
                access_schema=res_dict[resource]['schema']
                )



    # Shared data
    res.shared_data = [ './ala2.pdb','./simulate.py','./analyze.py' ]

    # Submit request for resources + wait till job becomes Active
    res.allocate(wait=True)

    try:

        # Run the first workload
        res.run(app)
    
    except Exception,ex:

        print 'Failed with error: ',ex

    finally:
        # Deallocate the resource
        res.deallocate()
