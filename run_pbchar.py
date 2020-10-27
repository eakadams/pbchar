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
    #updating table for new way of looking at values/errors
    dtype = [('beam','i4'),('mean_flux_ratio','f8'),
             ('median_flux_ratio','f8'),('std_flux_ratio','f8'),
             ('mean_error','f8'),('median_error','f8')]
    #two tables; one all vals, one 50%
    t = Table(data = np.zeros(N,dtype=dtype))
    t50 = Table(data = np.zeros(N,dtype=dtype))
    
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
        *ratio_vals = CB.get_ratio_vals()
        *ratio_vals_50 = CB.get_ratio_vals(level=0.5)
        
        #set trable row values
        t[bm] = (bm,ratio_vals)
        t50[bm] = (bm,ratio_vals_50)

    #write table out
    tablename = "ratio_values_{}.csv".format(args.matchfilebase)
    tablepath = os.path.join(args.path,tablename)
    tablename50 = "ratio_values_50_{}.csv".format(args.matchfilebase)
    tablepath50 = os.path.join(args.path,tablename)
    t.write(tablepath, format='csv',overwrite=True)
    t50.write(tablepath50, format='csv',overwrite=True)

    #make some plots
    fig, (ax1,ax2) = plt.subplots(1,2,figsize=(10,5))
    #plot ratio values
    ax1.scatter(t['bm'],t['mean_flux_ratio'],label='Mean flux ratio')
    ax1.scatter(t['bm'],t['median_flux_ratio'],label='Median flux ratio')
    ax1.scatter(t50['bm'],t50['mean_flux_ratio'],
                label='Mean flux ratio; >=50%')
    ax1.scatter(t50['bm'],t50['median_flux_ratio'],
                label='Median flux ratio; >=50%')
    ax1.set_xlabel('Beam')
    ax1.set_ylabel('Integrated flux ratio Apertif / NVSS')
    ax1.legend()

    #plot error/scatter values
    ax2.scatter(t['bm'],t['std_flux_ratio'],label='Scatter flux ratio')
    ax2.scatter(t['bm'],t['mean_error'],label='Mean error flux ratio')
    ax2.scatter(t['bm'],t['median_error'],label='Median error flux ratio')
    ax2.scatter(t50['bm'],t50['mean_error'],
                label='Mean error flux ratio; >=50%')
    ax2.scatter(t50['bm'],t50['median_error'],
                label='Median error flux ratio; >=50%')
    ax2.scatter(t50['bm'],t50['std_flux_ratio'],label='Scatter flux ratio; >=50%')
    ax2.set_xlabel('Beam')
    ax2.set_ylabel('Uncertainty Apertif/NVSS int flux')
    ax2.legend()

    #save figure!
    figpath = os.path.join(figdir,"ratio_vals.pdf")
    plt.savefig(figpath)
    plt.close('all')
    

