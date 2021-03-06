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
import matplotlib.pyplot as plt

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
parser.add_argument("--SN",default=None,type=float,
                    help="Min SN for filtering sources")
parser.add_argument("--mode",default="int",
                    help="'peak' or 'int', type of flux ratio")
parser.add_argument("--raterr",default=None,type=float,
                    help="Max error on flux ratio for filtering sources")
parser.add_argument("--size",default=None,type=float,
                     help="Max NVSS size for filtering sources")
parser.add_argument("--plots",default=False,type=bool,
                    help="Produce plots or not")
parser.add_argument("--output",default=None,
                    help="Additional string to add to output files")
args = parser.parse_args()


def beam_char(bm):
    """
    Make characterization plots for given beam
    """
    matchfile = os.path.join(args.path,"{0}_{1:02d}.csv".format(
        args.matchfilebase,bm))
    
    CB = pbchar.PB(matchfile,startdate=args.startdate,
                   enddate=args.enddate, SN=args.SN,
                   size = args.size, raterr=args.raterr,
                   mode = args.mode)
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

    #get a table with values for 2-D position plot
    dtype_pos = [('beam','i4'),('delta_ra','f8'),
                 ('delta_ra_cb','f8'), ('delta_dec','f8'),
                 ('delta_dec_cb','f8'), ('int_flux_ratio','f8')]
    t_pos = Table(data=np.zeros(N,dtype=dtype_pos))
    
    #run in serial per beam
    #weird issues with matplotlib / backend
    #and this is fast
    for bm in range(40):
        #get file
        matchfile = os.path.join(args.path,"{0}_{1:02d}.csv".format(
            args.matchfilebase,bm))
        #load PB object
        CB = pbchar.PB(matchfile,startdate=args.startdate,
                       enddate=args.enddate, SN=args.SN,
                       size = args.size, raterr=args.raterr,
                       mode = args.mode)
        #get plots
        #print(args.plots,type(args.plots))
        if args.plots:
            CB.go()
        #get scatter values for table
        (mean_ratio,median_ratio,std_ratio,
         mean_error,median_error) = CB.get_ratio_vals()
        (mean_ratio_50, median_ratio_50, std_ratio_50,
         mean_error_50, median_error_50) = CB.get_ratio_vals(level=0.5)
        
        #set trable row values
        t[bm] = (bm,mean_ratio,median_ratio,std_ratio,
                 mean_error,median_error)
        t50[bm] = (bm,mean_ratio_50, median_ratio_50,
                   std_ratio_50, mean_error_50, median_error_50)

    #get a limited table for docs and write that out
    t_docs = t['beam','median_flux_ratio','std_flux_ratio','median_error']
    t50_docs = t50['beam','median_flux_ratio','std_flux_ratio','median_error']
    tablename = "ratio_values_doc_{}.csv".format(args.matchfilebase)
    tablepath = os.path.join(args.path,tablename)
    t_docs.write(tablepath,format='csv',overwrite=True,
                 formats = {'median_flux_ratio': '4.2f',
                            'std_flux_ratio': '4.2f',
                            'median_error':'4.2f'})
    tablename50 = "ratio_values_50_doc_{}.csv".format(args.matchfilebase)
    tablepath50 = os.path.join(args.path,tablename50)
    t50_docs.write(tablepath50,format='csv',overwrite=True,
                   formats = {'median_flux_ratio': '4.2f',
                              'std_flux_ratio': '4.2f',
                              'median_error':'4.2f'})

    #write table out
    if args.output is None:
        tablename = "ratio_values_{}.csv".format(args.matchfilebase)
        tablename50 = "ratio_values_50_{}.csv".format(args.matchfilebase)
    else:
        tablename = "ratio_values_{0}_{1}.csv".format(args.matchfilebase,
                                                      args.output)
        tablename50 = "ratio_values_50_{0}_{1}.csv".format(args.matchfilebase,
                                                           args.output)
    tablepath = os.path.join(args.path,tablename)
    tablepath50 = os.path.join(args.path,tablename)
    t.write(tablepath, format='csv',overwrite=True)
    t50.write(tablepath50, format='csv',overwrite=True)

    #make some plots
    fig, (ax1,ax2) = plt.subplots(1,2,figsize=(10,5))
    #plot ratio values
    ax1.scatter(t['beam'],t['mean_flux_ratio'],
                label='Mean flux ratio: {:4.2f}'.format(np.median(
                    t['mean_flux_ratio'])))
    ax1.scatter(t['beam'],t['median_flux_ratio'],
                label='Median flux ratio: {:4.2f}'.format(np.median(
                    t['median_flux_ratio'])))
    ax1.scatter(t50['beam'],t50['mean_flux_ratio'],
                label='Mean flux ratio >=50% : {:4.2f}'.format(np.median(
                    t50['mean_flux_ratio'])))
    ax1.scatter(t50['beam'],t50['median_flux_ratio'],
                label='Median flux ratio >=50% : {:4.2f}'.format(np.median(
                    t50['median_flux_ratio'])))
    ax1.set_xlabel('Beam')
    ax1.set_ylabel('Integrated flux ratio Apertif / NVSS')
    ax1.legend()

    #plot error/scatter values
    ax2.scatter(t['beam'],t['std_flux_ratio'],
                label='Scatter flux ratio : {:4.2f}'.format(np.median(
                    t['std_flux_ratio'])))
    #ax2.scatter(t['beam'],t['mean_error'],
    #            label='Mean error flux ratio : {:4.2f}'.format(np.median(
    #                t['mean_error'])))
    ax2.scatter(t['beam'],t['median_error'],
                label='Median error flux ratio : {:4.2f}'.format(np.median(
                    t['median_error'])))
    ax2.scatter(t50['beam'],t50['std_flux_ratio'],
                label='Scatter flux ratio >=50% : {:4.2f}'.format(np.median(
                    t50['std_flux_ratio'])))
    #ax2.scatter(t50['beam'],t50['mean_error'],
    #            label='Mean error flux ratio >=50% : {:4.2f}'.format(np.median(
    #                t50['mean_error'])))
    ax2.scatter(t50['beam'],t50['median_error'],
                label='Median error flux ratio >=50% : {:4.2f}'.format(np.median(
                    t50['median_error'])))
    
    ax2.set_xlabel('Beam')
    ax2.set_ylabel('Uncertainty Apertif/NVSS int flux')
    ax2.legend()

    #save figure!
    if args.output is None:
        figpath = os.path.join(figdir,"ratio_vals_{}.pdf".format(args.matchfilebase))
    else:
        figpath = os.path.join(figdir,"ratio_vals_{0}_{1}.pdf".format(args.matchfilebase,
                                                                      args.output))
    plt.savefig(figpath)
    plt.close('all')
    

