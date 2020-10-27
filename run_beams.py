#python code to run beam in bulk

from __future__ import print_function

__author__  = "E.A.K. Adams"

"""
Run primary beam correction and 
NVSS cross-matching in bulk on a 
batch of beams.
Runs on all of DR1 for now
Can specify primary beam to use
but only 190912 and 191002 are implemented

Issues with processes dying during parallelization
Add try/except following
http://openkb.fr/Python_multiprocessing_pool_hangs_at_join

Future improvements:
- Options to run on other data
- Implement PB selection properly
- Add optional number of cores to argparse
"""

import beam as beam
import multiprocessing
import multiprocessing.pool
from astropy.io import ascii
from multiprocessing import Pool
import argparse
import os
import shutil

parser = argparse.ArgumentParser(description='Run primary beam correction and cross-matching')
parser.add_argument("--PB", help='Primary beam to use',
                    default='190912')
parser.add_argument("--ncores",help='Number of cores to use',
                    default=12, type=int)
parser.add_argument("--check",help="Check whether output exists before mapping/running",
                    default=True, type=bool)
parser.add_argument("--workingdir",help="Location to write tmp files",
                    default="/tank/adams/pbchar/tmp")
parser.add_argument("-outputdir",help="Location to write output",
                    default="/tank/adams/pbchar/")
args = parser.parse_args()

#get pbname and pbdir based on input
pbname = args.PB
if pbname == '190912':
    pbdir = '/tank/apertif/driftscans/DR1/190912/continuum'
elif pbname == '191002':
    pbdir='/tank/denes/drift_scan/191002/cont'
elif pbname == 'gpall':
    pbdir = '/data/kutkin/cbeams'
elif pbname == '200130':
    pbdir = '/tank/denes/drift_scan/200130/cont'
elif pbname == 'norm':
    pbdir = '/data/kutkin/cbeams'
elif pbname == 'nonorm':
    pbdir = '/data/kutkin/cbeams'
else:
    pbdir = '/tank/apertif/driftscans/DR1/190912/continuum'
    print('Using default DR PBs')


#cleanup working dir before starting
#avoids conflicts between runs
#can't run multiple at the same time!
if os.path.isdir(args.workingdir):
    #delete
    print("Deleting working directory before start")
    print(args.workingdir)
    shutil.rmtree(args.workingdir,ignore_errors = True)
    

#function to setup/run a beam
def run_beam(packed_args):
    #unpack  args
    bm,tid = packed_args
    #print("Initializing beam {0} for taskid {1}".format(bm,tid))
    B = beam.Beam(tid,bm,pbname=pbname,pbdir=pbdir,
                  masklevel=0.1,
                  outputdir = args.outputdir,
                  workingdir = args.workingdir)
    try:
        B.go()
    except Exception as e:
        print(("Processing failed for beam {0}, "
               "taskid {1}").format(bm,tid))
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
    pool = MyPool(args.ncores)
    #setup jobs
    jobs = []
    cont_obs = ascii.read('dr_year1_cont.csv')
    if args.check:
        for bm,tid in cont_obs['Beam','ObsID']:
            #first check if output exists
            beam_outputdir = os.path.join(args.outputdir,pbname,"{:02d}".format(bm))
            match_output = os.path.join(beam_outputdir,
                                        "{0}_matches.csv".format(tid))
            if not os.path.exists(match_output):
                jobs.append((bm,tid))
    else:
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
