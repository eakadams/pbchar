#python code to run beam in bulk

from __future__ import print_function

__author__  = "E.A.K. Adams"

"""
Run primary beam correction and 
NVSS cross-matching in bulk on a 
batch of beams.
Runs on all of DR1 for now
Can specify primary beam to use
but only default is implemented

Future improvements:
- Options to run on other data
- Parallelize
- Implement PB selection properly
- Add optional number of cores to argparse
"""

import beam as beam
from astropy.io import ascii
from multiprocessing import Pool
import argparse


parser = argparse.ArgumentParser(description='Run primary beam correction and cross-matching')
parser.add_argument("--PB", help='Primary beam to use',
                    default='190912')
args = parser.parse_args()

#get pbname and pbdir based on input
pbname = args.PB
if pbname == '190912':
    pbdir = '/tank/apertif/driftscans/DR1/190912/continuum'
elif pbname == '191002':
    pbdir='/tank/denes/drift_scan/191002/cont'
else:
    pbdir = '/tank/apertif/driftscans/DR1/190912/continuum'
    print('Using default DR PBs')
    

#function to setup/run a beam
def run_beam(packed_args):
    #unpack  args
    bm,tid = packed_args
    #print("Initializing beam {0} for taskid {1}".format(bm,tid))
    B = beam.Beam(tid,bm,pbname=pbname,pbdir=pbdir,
                  masklevel=0.1)
    B.go()
    #print("Ran beam {0} of taskid {1}".format(bm,tid))
    


#try to setup parallelization
if __name__ == '__main__':
    #get taskid/beam combination for DR1
    cont_obs = ascii.read('dr_year1_cont.csv')
    #test on just ten  entries
    #cont_obs = cont_obs[0:10]
    #set up pool
    #default 4 cores; will add as param
    pool  = Pool(6)
    #set up jobs list
    jobs = []
    for bm,tid in cont_obs['Beam','ObsID']:
        jobs.append((bm,tid))
    #print(jobs)
    #run jobs
    pool.map(run_beam,jobs)
    #close thing because i think I should
    pool.close()
    pool.join()
    
#for bm,taskid in cont_obs['Beam','ObsID']:
#    B = beam.Beam(taskid,bm,masklevel=0.1)
#    B.go()
