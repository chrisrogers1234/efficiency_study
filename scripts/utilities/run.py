#!/usr/env python

import copy
import subprocess
import time
import os
import shutil
import sys

class Run(object):
    def __init__(self, config):
        self.prefix = ['simulation']
        self.config = config
        self.jobs = config['jobs']
        self.n_procs = config['n_procs']
        self.extra_args = config['extra_args']
        self.delta_t = config['delta_t']
        self.max_t = config['max_t']
        self.logs = config['logs']
        self.processes = []
        self.make_log_dir()

    def make_log_dir(self):
        try:
            shutil.rmtree(self.logs)
        except OSError:
            pass
        os.makedirs(self.logs+"/sbatch")

    def run_many(self):
        print 'Running', len(self.jobs), 'jobs with cards', self.prefix,
        if self.is_scarf():
            print "on scarf"
        else:
            print "locally"
        print self.jobs
        timer = 0
        job_list = copy.deepcopy(self.jobs)

        self.processes = []
        while len(self.processes) > 0 or len(job_list) > 0:
            print '\r', round(timer/60., 1), '    ',
            self.poll()
            print len(self.processes),
            sys.stdout.flush()
            while len(self.processes) < self.n_procs and len(job_list) > 0:
                job = job_list.pop(0)
                for cards in self.prefix:
                    self.processes.append(self.run_one(job, cards))
            timer += self.delta_t
            if timer > self.max_t:
                print 'Out of time ... killing all jobs'
                self.kill_all()
                job_list = []
            time.sleep(self.delta_t)
        try:
            print subprocess.check_output(['squeue', '-u', 'scarf148'])
        except OSError:
            pass # not on scarf
        print "Done"

    def make_sbatch_file(self, unique_id, command):
        file_name = "logs/tmp/run_one_analysis_"+unique_id+".sbatch"
        fout = open(file_name, "w")
        print >> fout, "#!/bin/sh\n#SBATCH"
        for item in command:
            print >> fout, item,
        print >> fout
        return file_name

    def run_one(self, unique_id, cards):
        sbatch_name = cards.split(".")[0]
        sbatch_name = sbatch_name.split("/")[-1]
        unique_id = str(unique_id)
        log_name = self.logs+"/analysis_"+sbatch_name+"_"+unique_id+".log"
        command = ['python', 'scripts/bin/run_one_analysis.py', cards, unique_id]+self.extra_args
        if self.is_scarf():
            sbatch_filename = self.make_sbatch_file(unique_id, command)
            print "Command for sbatch", command
            command = ['sbatch',
                    '-n', '1',
                    '--time', '6-23', # minutes, 3 days 0 hrs
                    '-p', 'ibis',
                    '-o', log_name,
                    '-e', log_name,
                    sbatch_filename
                ]
            sbatch_name = self.logs+"/sbatch/"+sbatch_name+"_"+str(unique_id)+"_bsub.log"
            log_file = open(sbatch_name, 'w')
        else:
            log_file = open(log_name, 'w')
        subproc = subprocess.Popen(command, stdout=log_file, stderr=subprocess.STDOUT)
        print "Running command", command, "in process id", subproc.pid
        print 
        return subproc, sbatch_name

    @classmethod
    def is_scarf(cls):
        uname = subprocess.check_output(['uname', '-a'])
        return 'scarf.rl.ac.uk' in uname

    def get_bjob_number(self, bsub_name):
        line = open(bsub_name).readline()
        if line == "":
            return None
        bjob_number = line.split(' ')[3]
        return int(bjob_number)

    def poll(self, verbose = False):
        if verbose:
            print '\nPolling local'
        processes_update = self.poll_local(verbose)
        if self.is_scarf():
            if verbose:
                print 'Polling scarf'
            processes_update += self.poll_scarf(verbose)
        # remove duplicates
        processes_update = list(set(processes_update))
        self.processes = processes_update
            
    def poll_local(self, verbose):
        processes_update = []
        for proc, dir_name in self.processes:
            if proc.returncode == None:
                proc.poll()
                if proc.returncode == None:
                    processes_update.append((proc, dir_name))
            if verbose:
                print '   ', proc.pid, dir_name, proc.returncode
        return processes_update

    def poll_scarf(self, verbose):
        processes_update = []
        for proc, bsub_name in self.processes:
            bjob_number = self.get_bjob_number(bsub_name)
            output = subprocess.check_output(['squeue', '-u', 'scarf148'])
            okay = False
            for line in output.split('\n'):
                if str(bjob_number) in line:
                    processes_update.append((proc, bsub_name))
                    okay = True
                    break
            if not okay:
                print "Process with id", proc.pid, "and bsub log", bsub_name, "died"
            if verbose:
                print '   ', proc.pid, bjob_number
        return processes_update

    def kill_all(self):
        self.kill_all_local()
        if self.is_scarf():
            self.kill_all_scarf()
    
    def kill_all_local(self):
        for proc, dir_name in self.processes:
            if proc.returncode != None:
                continue # proc has returned
            pid_str = str(proc.pid)
            subprocess.check_output(['kill', '-9', pid_str])
            
    def kill_all_scarf(self):
        for proc, dir_name in self.processes:
            bjob_number = self.get_bjob_number(dir_name)
            try:
                output = subprocess.check_output(['scancel', str(bjob_number)])
            except Exception:
                pass

