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
            "name":         "mdrun",
            "description":  "MD simulation kernel",
            "arguments":   {"--tpr=":     
                        {
                            "mandatory": True,
                            "description": "Tpr file name"
                        },
                        "--rcon=":     
                        {
                            "mandatory": True,
                            "description": "Configuration file."
                        },
                        '--xtc=':
                        {
                            "mandatory": True,  
                            "description": "Trajectory file."
                        }
                    },
            "machine_configs": 
            {
                "*": {
                    "environment"   : None,
                    "pre_exec"      : None,
                    "executable"    : "mdrun",
                    "uses_mpi"      : False
                },
                "xsede.stampede":{
                    "environment"   : None,
                    "pre_exec"      : ['. /opt/apps/lmod/lmod/init/sh','module restore','module load boost','module load gromacs/5.1.2'],
                    "executable"    : "gmx mdrunp",
                    "uses_mpi"      : False
                },
                "local.localhost":{
                    "environment"   : None,
                    "pre_exec"      : ['export PATH=$PATH:/home/vivek91/modules/gromacs-5.1.3/build/bin'],
                    "executable"    : "gmx mdrun",
                    "uses_mpi"      : False
                },
                "xsede.comet":{
                    "environment"   : None,
                    "pre_exec"      : ['. /usr/share/Modules/init/sh','module load gromacs'],
                    "executable"    : "gmx mdrun",
                    "uses_mpi"      : False
                }
            }
    }


# ------------------------------------------------------------------------------
# 
class mdrun_kernel(KernelBase):

    # --------------------------------------------------------------------------
    #
    def __init__(self):
        """Le constructor.
        """
        super(mdrun_kernel, self).__init__(_KERNEL_INFO)


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
        arguments  = [  '-quiet', 
                        '-s', self.get_arg("--tpr="), 
                        '-cpi', 'state.cpt',
                        '-rcon', self.get_arg("--rcon="), 
                        '-x', self.get_arg("--xtc")
                    ]

        self._executable  = executable
        self._arguments   = arguments
        self._environment = cfg["environment"]
        self._uses_mpi    = cfg["uses_mpi"]
        self._pre_exec    = cfg["pre_exec"]

