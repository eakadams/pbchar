#python code to run pb char for all beams

from __future__ import print_function

__author__ = "E.A.K. Adams"

"""
This runs pbchar just for internal flux comparison
For a given PB name
"""

import pbchar as pbchar
import argparse
from multiprocessing import Pool
import os
import numpy as np
from astropy.table import Table, vstack
import matplotlib.pyplot as plt
from astropy.io import ascii
import matplotlib
import astropy.units as u
from scipy.stats import norm
import matplotlib.mlab as mlab

#get global level paths
this_dir,this_filename = os.path.split(__file__)
pbchardir = this_dir
filedir = os.path.join(pbchardir,"files")
figdir = os.path.join(pbchardir,"figures")

#set global params
colors = plt.rcParams['axes.prop_cycle'].by_key()['color']

def get_stats(array):
    """
    Get basic stats on array
    """
    median = np.nanmedian(array)
    mean = np.nanmedian(array)
    per16 = np.nanpercentile(array,16)
    per84 = np.nanpercentile(array,84)

    return [median, mean, per16, per84]

if __name__ == '__main__':

    #empty table list to hold internal comp tables
    tablelist = []
    
    #run per beam
    for bm in range(40):
        #get file
        matchfile = os.path.join('files',"gp_avg_orig_{0:02d}.csv".format(bm))
        #load PB object
        CB = pbchar.PB(matchfile)

        #get internal flux comparison tables
        comp_table = CB.internal_comp()
        tablelist.append(comp_table)
        


    #create a vstacked comparison table
    comp_table_all_beams = vstack(tablelist)
    print(comp_table_all_beams.colnames)
    print(len(comp_table_all_beams))

    #onlykeep SN>10
    #this is for all measurements
    #have I set SN to filter on Apertif and/or NVSS?
    #let's do both
    #first nvss
    nvss_sn = comp_table_all_beams['int_flux_nvss'] / comp_table_all_beams['int_flux_nvss_err']
    ind_nvss_sn = np.where(nvss_sn > 10)[0]
    comp_table_all_beams = comp_table_all_beams[ind_nvss_sn]
    #then apertif, both int and peak
    ap_int_sn = comp_table_all_beams['int_flux_ap'] / comp_table_all_beams['int_flux_ap_err']
    ind_ap_int_sn = np.where(ap_int_sn > 10)[0]
    comp_table_all_beams = comp_table_all_beams[ind_ap_int_sn]
    ap_peak_sn = comp_table_all_beams['peak_flux_ap'] / comp_table_all_beams['peak_flux_ap_err']
    ind_ap_peak_sn = np.where(ap_peak_sn > 10)[0]
    comp_table_all_beams = comp_table_all_beams[ind_ap_peak_sn]

    #and do size filtering. require NVSS majoraxis < 45"
    ind_size = np.where(comp_table_all_beams['maj_nvss'] < 45)[0]
    comp_table_all_beams = comp_table_all_beams[ind_size]

    #require at least three visits
    ind_3 = np.where(comp_table_all_beams['n_visits'] > 2)[0]
    comp_table_all_beams = comp_table_all_beams[ind_3]

    #after all filtering check length again
    print(len(comp_table_all_beams))

    #get normalized flux differences
    normdiff_int_flux_mean = ( ( comp_table_all_beams['int_flux_ap'] -
                                 comp_table_all_beams['mean_ap_int_flux'] ) /
                               comp_table_all_beams['mean_ap_int_flux'] )
    normdiff_peak_flux_mean = ( ( comp_table_all_beams['peak_flux_ap'] -
                                  comp_table_all_beams['mean_ap_peak_flux'] ) /
                                comp_table_all_beams['mean_ap_peak_flux'] )

    int_ratio = ( comp_table_all_beams['int_flux_ap'] /
                  comp_table_all_beams['mean_ap_int_flux'] )

    peak_ratio = ( comp_table_all_beams['peak_flux_ap'] /
                  comp_table_all_beams['mean_ap_peak_flux'] )

    #get diffs w/in 50% pb response
    ind_50 = np.where(comp_table_all_beams['pb_level'] >= 0.5)[0]
    normdiff_int_50 = normdiff_int_flux_mean[ind_50]
    normdiff_peak_50 = normdiff_peak_flux_mean[ind_50]
    int_ratio_50 = int_ratio[ind_50]
    peak_ratio_50 = int_ratio[ind_50]
    
    #make figures for internal flux comparison
    #just do a histogram

    fig, ((ax1) )= plt.subplots(1,1,figsize=(5,5))



    histbins = np.arange(-1.5,1.55,0.05)

    ax1.hist(int_ratio,
             bins = histbins,
             color=colors[0],
             histtype='step',
             label = 'Integrated flux')
    ax1.hist(int_ratio_50, bins=histbins,
             histtype='step',color=colors[0],
             linestyle='dashed',
             label = 'Integrated flux (>50% PB)')
    
    ax1.hist(peak_ratio,
             bins = histbins,
             color=colors[1],
             histtype='step',
             label = 'Peak flux')
    ax1.hist(peak_ratio_50, bins=histbins,
             histtype='step',color=colors[1],
             linestyle='dashed',
             label = 'Peak flux (>50% PB)')

    ax1.set_xlabel('Flux ratio (to mean)')
    

    ax1.set_yscale('log')
    ax1.set_ylabel('Counts')
    ax1.set_ylim([0.5,100000])

    ax1.legend(loc=2)
    
    stats_peak = get_stats(normdiff_peak_flux_mean)
    stats_peak_50 = get_stats(normdiff_peak_50)
    stats_int = get_stats(normdiff_int_flux_mean)
    stats_int_50 = get_stats(normdiff_int_50)

    stats_peak_ratio = get_stats(peak_ratio)
    stats_int_ratio = get_stats(int_ratio)
    stats_peak_ratio_50 = get_stats(peak_ratio_50)
    stats_int_ratio_50 = get_stats(int_ratio_50)

    print(('Integrated flux ratio (to the mean) is: '
           "{0:4.2f} (+{1:4.2f} - {2:4.2f})").format(
               stats_int_ratio[1],
               (stats_int_ratio[3] - stats_int_ratio[1]),
               (stats_int_ratio[1] - stats_int_ratio[2]) ) )
    
    print(('Integrated flux ratio (>50%) is: '
           "{0:4.2f} (+{1:4.2f} - {2:4.2f})").format(
               stats_int_ratio_50[1],
               (stats_int_ratio_50[3] - stats_int_ratio_50[1]),
               (stats_int_ratio_50[1] - stats_int_ratio_50[2]) ) )
    
    print(('Peak flux ratio (to the mean) is: '
           "{0:4.2f} (+{1:4.2f} - {2:4.2f})").format(
               stats_peak_ratio[1],
               (stats_peak_ratio[3] - stats_peak_ratio[1]),
               (stats_peak_ratio[1] - stats_peak_ratio[2]) ) )
    
    print(('Peak flux ratio (>50%) is: '
           "{0:4.2f} (+{1:4.2f} - {2:4.2f})").format(
               stats_peak_ratio_50[1],
               (stats_peak_ratio_50[3] - stats_peak_ratio_50[1]),
               (stats_peak_ratio_50[1] - stats_peak_ratio_50[2]) ) )



    figpath = os.path.join(figdir,"internal_flux_comp_hist.pdf")
    plt.savefig(figpath)
    plt.close()
