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

#get global level paths
this_dir,this_filename = os.path.split(__file__)
pbchardir = this_dir
filedir = os.path.join(pbchardir,"files")
figdir = os.path.join(pbchardir,"figures")

#get command line arguments
parser = argparse.ArgumentParser(description="Run characterization for all CBs")
parser.add_argument("matchfilebase", help="base name of csv file with w/ cross-matches") 
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
    #run in serial per beam
    #weird issues with matplotlib / backend
    #and this is fast
    for bm in range(40):
        beam_char(bm)
                    

    

