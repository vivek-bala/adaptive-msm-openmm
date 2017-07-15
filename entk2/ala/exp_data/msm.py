#!/usr/bin/env python

"""A kernel for OpenMM simulations
"""

__author__    = "Vivek <vivek.balasubramanian@rutgers.edu>"
__copyright__ = "Copyright 2017, http://radical.rutgers.edu"
__license__   = "MIT"

from copy import deepcopy

from radical.entk import NoKernelConfigurationError
from radical.entk import KernelBase

# ------------------------------------------------------------------------------
# 
_KERNEL_INFO = {
			"name":         "msm",
			"description":  "Execute openmm python scripte",
			"arguments":   {"--lag=":     
						{
							"mandatory": True,
							"description": "Lag time"
						},
					"--clusters=":
						{
							"mandatory": True,
							"description": "Number of clusters"
						},
                    "--pdb=":
                        {
                            "mandatory": True,
                            "description": "Pdb filename"
                        },
                    "--components=":
                        {
                            "mandatory": True,
                            "description": "Number of components"
                        },
                    "--stride=":
                        {
                            "mandatory": True,
                            "description": "Number of components"
                        },

					},
			"machine_configs": 
			{
				"xsede.stampede": {
					"environment"   : None,
					"pre_exec"      : ['export PATH=/home1/02734/vivek91/miniconda2/bin:$PATH','source activate msm_env'],
					"executable"    : "python",
					"uses_mpi"      : False
				},

                "xsede.supermic": {
                    "environment"   : None,
                    "pre_exec"      : ['export PATH=/home/vivek91/miniconda2/bin:$PATH','source activate openmm_env'],
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
		arguments  = ['analyze.py','--lag', self.get_arg('--lag='),
                        '--clusters', self.get_arg('--clusters='),
                        '--pdb', self.get_arg('--pdb='),
                        '--components', self.get_arg('--components='),
                        '--stride', self.get_arg('--stride=') ]

		self._executable  = executable
		self._arguments   = arguments
		self._environment = cfg["environment"]
		self._uses_mpi    = cfg["uses_mpi"]
		self._pre_exec    = cfg["pre_exec"]

