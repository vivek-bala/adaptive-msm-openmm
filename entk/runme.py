__author__    = "Vivek Balasubramanian <vivek.balasubramanian@rutgers.edu>"
__copyright__ = "Copyright 2016, http://radical.rutgers.edu"
__license__   = "MIT"

from radical.entk import EoP, AppManager, Kernel, ResourceHandle

from grompp import grompp_kernel
#from mdrun import mdrun_kernel

import argparse

ENSEMBLE_SIZE=4
PIPELINE_SIZE=1


def Test(EoP):

    def __init__(self, ensemble_size, pipeline_size):
        super(Test,self).__init__(ensemble_size, pipeline_size)

        def stage_1(self, instance):
            k1 = Kernel(name='grompp')
            k1.arguments = [  
                                "--mdp=grompp.mdp",
                                "--conf=equil{inst}.gro".format(instance),
                                "--top=topol.top",
                                "--out=topol.tpr"
                            ]        
            k1.cores=1

            k1.link_input_data = [
                                    '$SHARED/grompp.mdp',
                                    '$SHARED/equil{inst}.gro',
                                    '$SHARED/topol.top'
                                ]

            return k1


if __name__ == '__main__':

    # Create pattern object with desired ensemble size, pipeline size
    pipe = Test(ensemble_size=ENSEMBLE_SIZE, pipeline_size=PIPELINE_SIZE)

    # Create an application manager
    app = AppManager(name='MSM')

    # Register kernels to be used
    app.register_kernels(grompp_kernel)
    #app.register_kernels(mdrun_kernel)


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
                        '{0}/grompp.mdp'.format(path),
                        '{0}/topol.top'.format(path),
                    ]

    # Submit request for resources + wait till job becomes Active
    res.allocate(wait=True)

    # Run the given workload
    res.run(app)

    # Deallocate the resource
    res.deallocate()