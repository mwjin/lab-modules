"""
This script includes 'Job' class and functions for job schedulers (SGE and PBS).
* This script is a modified version of dandyrilla/baeklab/baeklab/base/jobsched.py.
* Original Author: Sukjun Kim
"""


class Job:
    """
    # TODO: make explanation
    """
    def __init__(self, **kwargs):
        self.id = None  # submission ID
        self.name = kwargs.get('name', 'Anonymous')
        self.queue = kwargs.get('queue', None)
        self.cmd = ''

        # attributes for job dependency
        self.prior_jobs = []  # a list of jobs that have higher priority to finish than this job


def qsub_sge():
    pass


def qsub_pbs():
    pass
