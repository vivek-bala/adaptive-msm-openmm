#!/usr/bin/env python

__author__    = "Vivek <vivek.balasubramanian@rutgers.edu>"
__copyright__ = "Copyright 2016, http://radical.rutgers.edu"
__license__   = "MIT"

from copy import deepcopy

from radical.entk import NoKernelConfigurationError
from radical.entk import KernelBase

# ------------------------------------------------------------------------------
# 
_KERNEL_INFO = {
            "name":         "msm",
            "description":  "MSM Analysis kernel",
            "arguments":   {"--micro=":     
                        {
                            "mandatory": True,
                            "description": "Number of microstates"
                        },
                        "--macro=":     
                        {
                            "mandatory": True,
                            "description": "Number of macrostates"
                        },
                        "--reference=":
                        {
                            "mandatory": True,
                            "description": "Reference pdb file"
                        },
                        "--grpname=":
                        {
                            "mandatory": True,
                            "description": "Group to study"
                        },
                        "--lag=":
                        {
                            "mandatory": True,
                            "description": "Lag time"
                        },
                        "--num_sims=":
                        {
                            "mandatory": True,
                            "description": "Number of simulations per macrostate"
                        }
#                        "--ensembles=":
#                        {
#                            "mandatory": True,
#                            "description": "Number of min simulations"
#                        }
                    },
            "machine_configs": 
            {
                "*": {
                    "environment"   : None,
                    "pre_exec"      : None,
                    "executable"    : "python",
                    "uses_mpi"      : False
                },
                "xsede.stampede":{
                    "environment"   : None,
                    "pre_exec"      : ['. /opt/apps/lmod/lmod/init/sh','module load python'],
                    "executable"    : "python",
                    "uses_mpi"      : False
                },
                "local.localhost":{
                    "environment"   : None,
                    "pre_exec"      : ['export PYTHONPATH=$PYTHONPATH:/usr/lib/python2.7/dist-packages:/home/vivek91/repos/adaptive-msm/tests:/home/vivek91/repos/adaptive-msm/tests/msmbuilder','export PATH=$PATH:/home/vivek91/modules/gromacs-5.1.3/build/bin:/usr/lib/python2.7/dist-packages'],
                    "executable"    : "python",
                    "uses_mpi"      : False
                },
                "xsede.comet":{
                    "environment"   : None,
                    "pre_exec"      : ['. /usr/share/Modules/init/sh','module load python'],
                    "executable"    : "python",
                    "uses_mpi"      : False
                }
            }
    }


# ------------------------------------------------------------------------------
# 
class msm_kernel(KernelBase):

    # --------------------------------------------------------------------------
    #
    def __init__(self):
        """Le constructor.
        """
        super(msm_kernel, self).__init__(_KERNEL_INFO)


    # --------------------------------------------------------------------------
    #
    def _bind_to_resource(self, resource_key):
        """(PRIVATE) Implements parent class method. 
        """
        if resource_key not in _KERNEL_INFO["machine_configs"]:
            if "*" in _KERNEL_INFO["machine_configs"]:
                # Fall-back to generic resource key
                resource_key = "*"
            else:
                raise NoKernelConfigurationError(kernel_name=_KERNEL_INFO["name"], resource_key=resource_key)

        cfg = _KERNEL_INFO["machine_configs"][resource_key]

        executable = cfg['executable']
        arguments  = [  
                        'MSMproject.py',
                        '--micro', self.get_arg("--micro="), 
                        '--macro', self.get_arg("--macro="), 
                        '--reference', self.get_arg("--reference="), 
                        '--grpname', self.get_arg("--grpname="), 
                        '--lag', self.get_arg("--lag="),
                        '--num_sims', self.get_arg("--num_sims=")
#                        '--ensembles', self.get_arg("--ensembles=")
                    ]

        self._executable  = executable
        self._arguments   = arguments
        self._environment = cfg["environment"]
        self._uses_mpi    = cfg["uses_mpi"]
        self._pre_exec    = cfg["pre_exec"]

