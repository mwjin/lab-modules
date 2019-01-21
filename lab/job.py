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
        self.prior_jobs = []  # a list of prior jobs in execution than this job
        self.hold_jid = kwargs.get('hold_jid', None)  # only used in SGE job scheduler
        self.priority = kwargs.get('priority', 0)  # a priority to execute (for dependency, only used in PBS scheduler)


def qsub_sge():
    pass


def qsub_pbs():
    pass
