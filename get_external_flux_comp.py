#python code to run pb char for all beams

from __future__ import print_function

__author__ = "E.A.K. Adams"

"""
This runs pbchar for external flux comparison
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
    #do the external flux comparison
    
    #empty table list to hold internal comp tables
    #do twice
    tablelist_orig = []
    tablelist_norm = []
    
    #run per beam
    for bm in range(40):
        #get file
        orig = os.path.join('files',"gp_avg_orig_{:02d}.csv".format(bm))
        norm = os.path.join('files',"gp_avg_norm_{:02d}.csv".format(bm))

        #load PB object
        #remove filtering from there, and do manually here if I want/need
        CBnorm = pbchar.PB(norm)
        CBorig = pbchar.PB(orig)
        
        #get external flux comparison tables
        orig_match = CBorig.matches
        orig_match['beam'] = np.full(len(orig_match),bm)
        norm_match = CBnorm.matches
        norm_match['beam'] = np.full(len(norm_match),bm)
    
        tablelist_orig.append(orig_match)
        tablelist_norm.append(norm_match)
    


    #create a vstacked comparison table
    orig_all = vstack(orig_match)
    norm_all = vstack(norm_match)

    #only keep SN > 10
    #this is for all measurements
    print(len(orig_all),len(norm_all))
    nvss_sn_orig = orig_all['int_flux_nvss'] / orig_all['int_flux_nvss_err']
    ind_nvss_orig = np.where(nvss_sn_orig > 10)[0]
    orig_all = orig_all[ind_nvss_orig]

    nvss_sn_norm = norm_all['int_flux_nvss'] / norm_all['int_flux_nvss_err']
    ind_nvss_norm = np.where(nvss_sn_norm > 10)[0]
    norm_all = norm_all[ind_nvss_norm]

    int_sn_orig = orig_all['int_flux_ap'] / orig_all['int_flux_ap_err']
    ind_int_orig = np.where(int_sn_orig > 10)[0]
    orig_all = orig_all[ind_int_orig]

    int_sn_norm = norm_all['int_flux_ap'] / norm_all['int_flux_ap_err']
    ind_int_norm = np.where(int_sn_norm > 10)[0]
    norm_all = norm_all[ind_int_norm]

    peak_sn_orig = orig_all['peak_flux_ap'] / orig_all['peak_flux_ap_err']
    ind_peak_orig = np.where(peak_sn_orig > 10)[0]
    orig_all = orig_all[ind_peak_orig]

    peak_sn_norm = norm_all['peak_flux_ap'] / norm_all['peak_flux_ap_err']
    ind_peak_norm = np.where(peak_sn_norm > 10)[0]
    norm_all = norm_all[ind_peak_norm]

    #also only keep if NVSS majoraxis < 45"
    ind_size_orig = np.where(orig_all['maj_nvss'] < 45)[0]
    ind_size_norm = np.where(norm_all['maj_nvss'] < 45)[0]

    norm_all = norm_all[ind_size_norm]
    orig_all = orig_all[ind_size_orig]

    print(len(orig_all),len(norm_all))

    #get nvss flux ratios
    int_ratio_orig = orig_all['int_flux_ap'] / orig_all['int_flux_nvss']
    peak_ratio_orig = orig_all['peak_flux_ap'] / orig_all['int_flux_nvss']

    int_ratio_norm = norm_all['int_flux_ap'] / norm_all['int_flux_nvss']
    peak_ratio_norm = norm_all['peak_flux_ap'] / norm_all['int_flux_nvss']

    print(orig_all.colnames)

   
    #get flux ratios at >=50% of pb response
    ind_orig_50 = np.where(orig_all['pb_level'] >= 0.5)[0]
    ind_norm_50 = np.where(norm_all['pb_level'] >= 0.5)[0]

    int_ratio_orig_50 = int_ratio_orig[ind_orig_50]
    int_ratio_norm_50 = int_ratio_norm[ind_norm_50]
    peak_ratio_orig_50 = peak_ratio_orig[ind_orig_50]
    peak_ratio_norm_50 = peak_ratio_norm[ind_norm_50]


    print(len(int_ratio_norm),len(int_ratio_norm_50),len(peak_ratio_norm),
          len(peak_ratio_norm_50))
    
    #make histogram of flux ratios

    fig, ((ax1, ax2) )= plt.subplots(1,2,figsize=(10,5) )

    histbins = np.arange(0.3,1.95,0.05)

    
    #get stats for labels
    stats_peak_orig = get_stats(peak_ratio_orig)
    stats_int_orig = get_stats(int_ratio_orig)
    stats_peak_orig_50 = get_stats(peak_ratio_orig_50)
    stats_int_orig_50 = get_stats(int_ratio_orig_50)


    
    ax1.hist(int_ratio_orig, bins=histbins,
             histtype='step', color=colors[0],
             linestyle='solid',
             label = ('Integrated flux ratio')
             )
    
    ax1.hist(int_ratio_orig_50, bins=histbins,
             histtype='step', color=colors[0],
             linestyle='dashed',
             label = ('Integrated flux ratio (>50% PB)')
    )

    
    
    ax1.hist(peak_ratio_orig, bins=histbins,
             histtype='step', color=colors[1],
             linestyle='solid',
             label = ('Peak flux ratio')
    )
    

    
    ax1.hist(peak_ratio_orig_50, bins=histbins,
             histtype='step', color=colors[1],
             linestyle='dashed',
             label = ('Peak flux ratio (>50% PB)')
             )

    ax1.set_xlabel('Apertif / NVSS integrated flux ratio, Original GPR')
    

    ax1.set_yscale('log')
    ax1.set_ylabel('Counts')

    ax1.set_ylim([0.5,5000])

    ax1.legend()


    print('For original GPR beams:')

    print(('Mean integrated flux ratio is: '
           "{0:4.2f} (+{1:4.2f} - {2:4.2f})").format(
               stats_int_orig[1],
               (stats_int_orig[3] - stats_int_orig[1]),
               (stats_int_orig[1] - stats_int_orig[2]) ) )
    
    print(('Mean integrated flux ratio (>50%) is: '
           "{0:4.2f} (+{1:4.2f} - {2:4.2f})").format(
               stats_int_orig_50[1],
               (stats_int_orig_50[3] - stats_int_orig_50[1]),
               (stats_int_orig_50[1] - stats_int_orig_50[2]) ) )

    print(('Mean peak flux ratio is: '
           "{0:4.2f} (+{1:4.2f} - {2:4.2f})").format(
               stats_peak_orig[1],
               (stats_peak_orig[3] - stats_peak_orig[1]),
               (stats_peak_orig[1] - stats_peak_orig[2]) ) )
    
    print(('Mean peak flux ratio (>50%) is: '
           "{0:4.2f} (+{1:4.2f} - {2:4.2f})").format(
               stats_peak_orig_50[1],
               (stats_peak_orig_50[3] - stats_peak_orig_50[1]),
               (stats_peak_orig_50[1] - stats_peak_orig_50[2]) ) )
    
  

    #now repeat for normalized beams, ax2

    #get stats for labels
    stats_peak_norm = get_stats(peak_ratio_norm)
    stats_int_norm = get_stats(int_ratio_norm)
    stats_peak_norm_50 = get_stats(peak_ratio_norm_50)
    stats_int_norm_50 = get_stats(int_ratio_norm_50)


    
    ax2.hist(int_ratio_norm, bins=histbins,
             histtype='step', color=colors[0],
             linestyle='solid',
             label = ('Integrated flux ratio')
             )
    
    ax2.hist(int_ratio_norm_50, bins=histbins,
             histtype='step', color=colors[0],
             linestyle='dashed',
             label = ('Integrated flux ratio (>50% PB)')
    )

    
    
    ax2.hist(peak_ratio_norm, bins=histbins,
             histtype='step', color=colors[1],
             linestyle='solid',
             label = ('Peak flux ratio')
    )
    

    
    ax2.hist(peak_ratio_norm_50, bins=histbins,
             histtype='step', color=colors[1],
             linestyle='dashed',
             label = ('Peak flux ratio (>50% PB)')
             )

    ax2.set_xlabel('Apertif / NVSS integrated flux ratio, Normalized GPR')
    

    ax2.set_yscale('log')
    ax2.set_ylabel('Counts')

    ax2.set_ylim([0.5,5000])

    ax2.legend()


    print('For normalized GPR beams:')

    print(('Mean integrated flux ratio is: '
           "{0:4.2f} (+{1:4.2f} - {2:4.2f})").format(
               stats_int_norm[1],
               (stats_int_norm[3] - stats_int_norm[1]),
               (stats_int_norm[1] - stats_int_norm[2]) ) )
    
    print(('Mean integrated flux ratio (>50%) is: '
           "{0:4.2f} (+{1:4.2f} - {2:4.2f})").format(
               stats_int_norm_50[1],
               (stats_int_norm_50[3] - stats_int_norm_50[1]),
               (stats_int_norm_50[1] - stats_int_norm_50[2]) ) )

    print(('Mean peak flux ratio is: '
           "{0:4.2f} (+{1:4.2f} - {2:4.2f})").format(
               stats_peak_norm[1],
               (stats_peak_norm[3] - stats_peak_norm[1]),
               (stats_peak_norm[1] - stats_peak_norm[2]) ) )
    
    print(('Mean peak flux ratio (>50%) is: '
           "{0:4.2f} (+{1:4.2f} - {2:4.2f})").format(
               stats_peak_norm_50[1],
               (stats_peak_norm_50[3] - stats_peak_norm_50[1]),
               (stats_peak_norm_50[1] - stats_peak_norm_50[2]) ) )

    print(('Standard deviation of '
           'integrated flux ratio is '
           '{0:4.2f}').format(np.std(int_ratio_norm)))
    
    
    figpath = os.path.join(figdir,"external_flux_comp_hist.pdf")

    plt.savefig(figpath)
    plt.close()

    
