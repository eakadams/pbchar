#python code to run beam in bulk

from __future__ import print_function

import beam as beam
from astropy.io import ascii
from multiprocessing import Pool

#get taskid/beam combination
#read in data to be released
cont_obs = ascii.read('dr_year1_cont.csv')

#test on just ten  entries
#cont_obs = cont_obs[0:10]

#try to setup parallelization
#in order to do this
#skipping parallelization for now
#want to figure out / learn how to do
#but I think I can brute force for tonight
#and do this when fresh
for bm,taskid in cont_obs['Beam','ObsID']:
    B = beam.Beam(taskid,bm,masklevel=0.1)
    B.go()
