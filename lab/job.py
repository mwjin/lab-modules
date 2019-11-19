"""
This script includes 'Job' class and functions for executing job schedulers (SGE, PBS, and sbatch).
Currently only one job per core is possible.
* This script is a modified version of dandyrilla/baeklab/baeklab/base/jobsched.py.
* Original Author: Sukjun Kim
"""

import os
import subprocess


class Job:
    """
    A class that represents one job
    """
    def __init__(self, name, cmd, **kwargs):
        """
        :param name: a name of an object of the 'Job' class
        :param cmd: a command line of this job
        :param kwargs:
            - hold_jid: a regular expression of prior job names (only for SGE job scheduler)
            - priority: an integer that represent a priority of this job for job dependency (only for PBS job scheduler)
        """
        self.name = name
        self.cmd = cmd

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
        echo_proc = subprocess.Popen(('echo', job.cmd), stdout=subprocess.PIPE)

        if job.hold_jid is None:
            qsub_proc = subprocess.Popen(('qsub', '-cwd', '-j', 'y', '-o', log_path, '-q', queue, '-N', job.name),
                                         stdin=echo_proc.stdout, stdout=subprocess.DEVNULL)
        else:
            qsub_proc = subprocess.Popen(('qsub', '-cwd', '-j', 'y', '-o', log_path, '-q', queue, '-N', job.name,
                                          '-hold_jid', job.hold_jid), stdin=echo_proc.stdout, stdout=subprocess.DEVNULL)

        echo_proc.wait()
        qsub_proc.wait()

        print('[LOG] Job (PID:%s) \'%s\' is submitted to the queue \'%s\'.' % (qsub_proc.pid, job.name, queue))


def qsub_pbs(jobs, queue, log_dir):
    """
    Throw jobs using PBS job scheduler
    :param jobs: a list of 'Job' objects
    :param queue: a queue name
    :param log_dir: a directory of logs for the jobs
    """
    os.makedirs(log_dir, exist_ok=True)
    priority_to_jobs = {}  # key: an integer that represents a priority of a job, value: a list of jobs

    for job in jobs:
        if job.priority not in priority_to_jobs:
            priority_to_jobs[job.priority] = []

        priority_to_jobs[job.priority].append(job)

    priorities = list(priority_to_jobs.keys())
    priorities.sort()

    prior_job_ids = []  # job IDs that have a higher priority

    for priority in priorities:
        curr_jobs = priority_to_jobs[priority]
        prior_job_ids_str = ':'.join(prior_job_ids)
        depend_opt = 'depend=afterok:%s' % prior_job_ids_str
        prior_job_ids = []  # reset

        for job in curr_jobs:
            log_path = '%s/%s.txt' % (log_dir, job.name)

            echo_proc = subprocess.Popen(('echo', job.cmd), stdout=subprocess.PIPE)
            qsub_proc = subprocess.Popen(('qsub', '-j', 'oe', '-o',  log_path, '-q', queue, '-N', job.name,
                                          '-W', depend_opt), stdin=echo_proc.stdout, stdout=subprocess.DEVNULL)

            echo_proc.wait()
            qsub_proc.wait()

            job_id = qsub_proc.pid
            prior_job_ids.append(job_id)

            print('[LOG] Job (PID:%s) \'%s\' is submitted to the queue \'%s\'.' % (job_id, job.name, queue))


def sbatch(jobs, partition, log_dir):
    """
    Throw jobs using SLURM job scheduler
    Currently setting the dependency between jobs is not available.
    :param jobs: a list of 'Job' objects
    :param partition: a partition name (equal with the 'queue; in SGE)
    :param log_dir: a directory of logs for the jobs
    """
    os.makedirs(log_dir, exist_ok=True)

    for job in jobs:
        log_path = '%s/%s.txt' % (log_dir, job.name)
        qsub_proc = subprocess.Popen(['sbatch', '-o', log_path, '-e', log_path, '-p', partition,
                                      '-J', job.name, f'--wrap={job.cmd}'], stdout=subprocess.DEVNULL)
        qsub_proc.wait()

        print('[LOG] Job (PID:%s) \'%s\' is submitted to the partition \'%s\'.' % (qsub_proc.pid, job.name, partition))

