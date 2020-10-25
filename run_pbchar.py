#python code to run pb char for all beams

from __future__ import print_function

__author__ = "E.A.K. Adams"

"""
Run pbchar code for all compound beams
For a given PB name
"""

import pbchar as pbchar
import argparse
from multiprocessing import Pool
import os
import numpy as np
from astropy.table import Table

#get global level paths
this_dir,this_filename = os.path.split(__file__)
pbchardir = this_dir
filedir = os.path.join(pbchardir,"files")
figdir = os.path.join(pbchardir,"figures")

#get command line arguments
parser = argparse.ArgumentParser(description="Run characterization for all CBs")
parser.add_argument("matchfilebase", help=("base name of csv file with w/ cross-matches, "
                                           "e..g, 'gpall'")) 
parser.add_argument("--startdate",help="Start date for plotting; YYMMDD",
                    default=None)
parser.add_argument("--enddate",help="End date for plotting; YYMMDD",
                    default=None)
parser.add_argument("--path",help=("Path to cross-matchfiles; "
                                   "default 'pbchar/files'"),
                    default=filedir)
args = parser.parse_args()


def beam_char(bm):
    """
    Make characterization plots for given beam
    """
    matchfile = os.path.join(args.path,"{0}_{1:02d}.csv".format(
        args.matchfilebase,bm))
    
    CB = pbchar.PB(matchfile,startdate=args.startdate,
                   enddate=args.enddate)
    CB.go()


if __name__ == '__main__':
    #setup up table to hold scatter output
    #40 rows, one for each beam
    N = 40
    #set up columns
    dtype = [('beam','i4'),('mean_flux_ratio','f8'),
             ('sc_N_50','f8'),('sc_N_mean','f8'),
             ('sc_bin_50','f8'),('sc_bin_mean','f8')]
    t = Table(data = np.zeros(N,dtype=dtype))
    
    #run in serial per beam
    #weird issues with matplotlib / backend
    #and this is fast
    for bm in range(40):
        #get file
        matchfile = os.path.join(args.path,"{0}_{1:02d}.csv".format(
            args.matchfilebase,bm))
        #load PB object
        CB = pbchar.PB(matchfile,startdate=args.startdate,
                       enddate=args.enddate)
        #get plots
        CB.go()
        #get scatter values for table
        (sc_N_50, sc_bin_50, sc_N_mean,
         sc_bin_mean) = CB.get_scatter_pblev(0.5,xbins = np.arange(0.15,1.1,0.1),
                                             mode = 'int')
        #get mean ratio value for table
        mean_flux_ratio = CB.get_mean_ratio(mode='int')
        
        #set trable row values
        t[bm] = (bm,mean_flux_ratio,
                 sc_N_50,sc_N_mean,
                 sc_bin_50,sc_bin_mean)

    #write table out
    tablename = "flux_ratio_mean_scatter_{}.csv".format(args.matchfilebase)
    tablepath = os.path.join(args.path,tablename)
    t.write(tablepath, format='csv',overwrite=True,
            formats = {'mean_flux_ratio': '4.2f',
                       'sc_N_50':'4.2f','sc_N_mean':'4.2f',
                       'sc_bin_50':'4.2f','sc_bin_mean':'4.2f'})
        
    
                    

    

