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
- Implement PB selection properly
- Add optional number of cores to argparse
"""

import beam as beam
import multiprocessing
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
#pyBDSF also uses pool so have to try something
#to deal with that:
#https://stackoverflow.com/questions/6974695/python-process-pool-non-daemonic
class NoDaemonProcess(multiprocessing.Process):
    #make 'daemon' attribute always return False
    def _get_daemon(self):
        return False
    def _set_daemon(self,value):
        pass
    daemon = property(_get_daemon, _set_daemon)

class MyPool(multiprocessing.pool.Pool):
    Process = NoDaemonProcess

def work():
    """
    Put work here, 
    pullint out of main test
    """
    #setup pool
    pool = MyPool(6)
    #setup jobs
    jobs = []
    for bm,tid in cont_obs['Beam','ObsID']:
        jobs.append((bm,tid))
    #run pool
    pool.map(run_beam,jobs)
    #close and join
    pool.close()
    pool.join()
    
    
if __name__ == '__main__':
   #do  work
    work()
