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
            "name":         "traj_collect",
            "description":  "Trajectory collection and processing to produce data required for MSMProject",
            "arguments":   {"--xtc=":     
                        {
                            "mandatory": True,
                            "description": "Trajectory files"
                        },
                        "--xtc_nopbc=":     
                        {
                            "mandatory": True,
                            "description": "Trajectory files without PBC"
                        },
                        "--system=":     
                        {
                            "mandatory": True,
                            "description": "System being monitored"
                        },
                        "--reference=":
                        {
                            "mandatory": True,
                            "description": "Reference filename"
                        },
                        "--lh5=":
                        {
                            "mandatory": True,
                            "description": "lh5 filename"
                        },
                        "--tpr=":
                        {
                            "mandatory": True,
                            "description": "tpr filename"
                        }
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
                    "pre_exec"      : ['export PYTHONPATH=$PYTHONPATH:/usr/lib/python2.7/dist-packages:/home/vivek91/repos/adaptive-msm/tests:/home/vivek91/repos/adaptive-msm/tests/msmbuilder','export PATH=$PATH:/home/vivek91/modules/gromacs-5.1.3/build/bin/usr/lib/python2.7/dist-packages:'],
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
class traj_collect_kernel(KernelBase):

    # --------------------------------------------------------------------------
    #
    def __init__(self):
        """Le constructor.
        """
        super(traj_collect_kernel, self).__init__(_KERNEL_INFO)


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
                        'pre_analysis.py',
                        '--xtc', self.get_arg("--xtc="),
                        '--system', self.get_arg("--system="), 
                        '--xtc_nopbc', self.get_arg("--xtc_nopbc="),
                        '--reference', self.get_arg("--reference="),
                        '--lh5', self.get_arg("--lh5="),
                        '--tpr', self.get_arg("--tpr=")
                    ]

        self._executable  = executable
        self._arguments   = arguments
        self._environment = cfg["environment"]
        self._uses_mpi    = cfg["uses_mpi"]
        self._pre_exec    = cfg["pre_exec"]

