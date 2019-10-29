#!/usr/bin/env python
#
# fast_ep_shelxe ->
# 
# code to run shelxe and manage the jobs - this will be called from within
# a multiprocess task so life is easier if the input is provided in the form
# of a dictionary with the information we need.

import os
import sys
import time
import multiprocessing
import shutil

os.path.dirname(os.path.abspath(__file__))
fast_ep_lib = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                           '../lib/')

if not fast_ep_lib in sys.path:
    sys.path.append(fast_ep_lib)

from run_job import run_job, run_job_cluster, is_cluster_job_finished

def run_shelxe_cluster(_settings):
    '''Run shelxe on cluster with settings given in dictionary, containing:

    nsite - number of sites
    solv - solvent fraction
    hand - original or inverted
    wd - working directory'''

    nsite = _settings['nsite']
    solv = _settings['solv']
    hand = _settings['hand']
    wd = _settings['wd']

    if hand == 'original':
        job_id = run_job_cluster(
#            'shelxe', ['sad', 'sad_fa', '-h%d' % nsite, '-z',
            'shelxe', ['sad', 'sad_fa', '-h%d' % nsite,
                       '-s%f' % solv, '-m20'], [], wd, 1, timeout = 600)
    else:
        job_id = run_job_cluster(
#            'shelxe', ['sad', 'sad_fa', '-h%d' % nsite, '-z',
            'shelxe', ['sad', 'sad_fa', '-h%d' % nsite,
                       '-s%f' % solv, '-m20', '-i'], [], wd, 1, timeout = 600)
        
    while not is_cluster_job_finished(job_id):
        time.sleep(1)

    return

def run_shelxe_local(_settings):
    '''Run shelxe locally with settings given in dictionary, containing:

    nsite - number of sites
    solv - solvent fraction
    hand - original or inverted
    wd - working directory'''

    nsite = _settings['nsite']
    solv = _settings['solv']
    hand = _settings['hand']
    wd = _settings['wd']

    if hand == 'original':
#        job_output = run_job('shelxe', ['sad', 'sad_fa', '-h%d' % nsite, '-z',
        job_output = run_job('shelxe', ['sad', 'sad_fa', '-h%d' % nsite,
                                        '-s%f' % solv, '-m20'], [], wd)
    else:
#        job_output = run_job('shelxe', ['sad', 'sad_fa', '-h%d' % nsite, '-z',
        job_output = run_job('shelxe', ['sad', 'sad_fa', '-h%d' % nsite,
                                        '-s%f' % solv, '-m20', '-i'], [], wd)

    return

def run_shelxe_multiprocess(_settings):
    '''Run multiple shelxe processes with various solvent contents and
    selection of hands'''

    wd_base = _settings['wd']
    nsite_real = _settings['nsite']

    solvent_fractions = [0.25 + 0.05 * j for j in range(11)]
    jobs = []
    
    for solvent_fraction in solvent_fractions:
        wd = os.path.join(wd_base, '%.2f' % solvent_fraction)
        if not os.path.exists(wd):
            os.makedirs(wd)
        shutil.copyfile(os.path.join(wd, '../../../../sad.hkl'),
                        os.path.join(wd, 'sad.hkl'))
        for ending in 'lst', 'pdb', 'res', 'hkl':
            shutil.copyfile(os.path.join(wd, '../sad_fa.%s' % ending),
                        os.path.join(wd, 'sad_fa.%s' % ending))

        jobs.append({'nsite':nsite_real, 'solv':solvent_fraction,
                     'hand':'original', 'wd':wd})
        jobs.append({'nsite':nsite_real, 'solv':solvent_fraction,
                     'hand':'inverted', 'wd':wd})
    
    p = multiprocessing.Pool(22)
    p.map(run_shelxe_local, jobs)

if __name__ == '__main__':
    dirbase = sys.argv[1]
    nsite = int(sys.argv[2])
    run_shelxe_multiprocess({'wd':dirbase,'nsite':nsite})
