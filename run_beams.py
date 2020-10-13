#python code to run beam in bulk

from __future__ import print_function

import beam as beam
from astropy.io import ascii

#get taskid/beam combination
#read in data to be released
cont_obs = ascii.read('dr_year1_cont.csv')

#to start with, only do first ten as test
cont_obs = cont_obs[0:10]

for beam,taskid in cont_obs['Beam','ObsID']:
    B = beam.Beam(taskid,beam,masklevel=0.1)
    B.go()
