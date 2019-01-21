"""
This script includes 'Job' class and functions for job schedulers (SGE and PBS).
* This script is a modified version of dandyrilla/baeklab/baeklab/base/jobsched.py.
* Original Author: Sukjun Kim
"""

import os
import sys
import subprocess


class Job:
    """
    A class that represents one job
    """
    def __init__(self, name, **kwargs):
        self.name = name
        self.cmd = kwargs.get('cmd', '')

        # attributes for job dependency
        self.hold_jid = kwargs.get('hold_jid', None)  # regex for names of prior jobs (only used in SGE job scheduler)
        self.priority = kwargs.get('priority', 0)  # a priority to execute (for dependency, only used in PBS scheduler)


def qsub_sge(jobs, queue, log_dir):
    """
    Throw jobs using Sun Grid Engine (SGE) job scheduler
    :param jobs: a list of 'Job' objects
    :param queue: a queue name
    :param log_dir: a directory of logs for the jobs
    """
    os.makedirs(log_dir, exist_ok=True)

    for job in jobs:
        log_path = '%s/%s.txt' % (log_dir, job.name)
        qsub_cmd = 'qsub -cwd -j y -o %s -q %s -N %s' % (log_path, queue, job.name)

        if job.hold_jid is not None:
            qsub_cmd += ' -hold_jid %s' % job.hold_jid

        echo_proc = subprocess.Popen(('echo', job.cmd), stdout=subprocess.PIPE)
        subprocess.Popen(('qsub', '-cwd', '-j', 'y', '-o', log_path, '-q', queue, '-N', job.name),
                         stdin=echo_proc.stdout, stdout=sys.stdout, stderr=sys.stderr)
        echo_proc.wait()


def qsub_pbs(jobs, queue, log_dir):
    """
    Throw jobs using PBS job scheduler
    :param jobs: a list of 'Job' objects
    :param queue: a queue name
    :param log_dir: a directory of logs for the jobs
    """
    priority_to_jobs = {}  # key: an integer that represents a priority of a job, value: a list of jobs

    for job in jobs:
        if job.priority not in priority_to_jobs:
            priority_to_jobs[job.priority] = []

        priority_to_jobs[job.priority].append(job)

    priorities = list(priority_to_jobs.keys())

    if len(priorities) == 1:  # no dependency
        pass
    else:
        priorities.sort()
        prior_job_ids = []  # job IDs that have a higher priority

        for priority in priorities:
            curr_jobs = priority_to_jobs[priority]
            prior_job_ids_str = ':'.join(prior_job_ids)
            depend_opt = 'depend=afterok:{%s}' % prior_job_ids_str
            prior_job_ids = []  # reset

            for job in curr_jobs:
                log_path = '%s/%s.txt' % (log_dir, job.name)

                echo_proc = subprocess.Popen(('echo', job.cmd), stdout=subprocess.PIPE)
                qsub_proc = subprocess.Popen(('qsub', '-j', 'oe', '-o',  log_path, '-q', queue, '-N', job.name,
                                              '-W', depend_opt), stdin=echo_proc.stdout)
                echo_proc.wait()
                job_id = qsub_proc.pid
                prior_job_ids.append(job_id)

                print('[LOG] Job %s \'%s\' is submitted to the queue \'%s\'.' % (job_id, job.name, queue))
