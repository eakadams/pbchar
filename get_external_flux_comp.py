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
    peak_ratio_orig_50 = peak_ratio_orig[ind_orig_50]

    int_ratio_norm_50 = int_ratio_norm[ind_norm_50]
    peak_ratio_norm_50 = peak_ratio_norm[ind_norm_50]

    #make histogram of flux ratios

    fig, ((ax1, ax2) )= plt.subplots(1,2,figsize=(10,5) )

    #get stats for labels
    stats_peak_orig = get_stats(peak_ratio_orig)
    stats_int_orig = get_stats(int_ratio_orig)
    stats_peak_orig_50 = get_stats(peak_ratio_orig_50)
    stats_int_orig_50 = get_stats(int_ratio_orig_50)

    labels_orig = [ ('Peak flux; mean = {0:4.2f} '
                     '(+{1:4.2f} -{2:4.2f})').format(stats_peak_orig[1],
                                                     stats_peak_orig[2],
                                                     stats_peak_orig[3]),
                    ('Int flux; mean = {0:4.2f} '
                     '(+{1:4.2f} -{2:4.2f})').format(stats_int_orig[1],
                                                     stats_int_orig[2],
                                                     stats_int_orig[3]),
                    ('Peak flux (>50% PB); mean = {0:4.2f} '
                     '(+{1:4.2f} -{2:4.2f})').format(stats_peak_orig_50[1],
                                                     stats_peak_orig_50[2],
                                                     stats_peak_orig_50[3]),
                    ('Int flux (>50% PB); mean = {0:4.2f} '
                     '(+{1:4.2f} -{2:4.2f})').format(stats_int_orig_50[1],
                                                     stats_int_orig_50[2],
                                                     stats_int_orig_50[3])
    ]

    histbins = np.arange(0.3,1.95,0.05)
    
    ax1.hist(int_ratio_orig, bins=histbins,
             histtype='step', color=colors[0],
             linestyle='solid',
             label = ('Int flux; mean = {0:4.2f} '
                     '(+{1:4.2f} - '
                      '{2:4.2f})').format(stats_int_orig[1],
                                          (stats_int_orig[3]-
                                           stats_int_orig[1]),
                                          (stats_int_orig[1] -
                                           stats_int_orig[2]))
             )
    
    ax1.hist(int_ratio_orig_50, bins=histbins,
             histtype='step', color=colors[0],
             linestyle='dashed',
             label = ('Int flux (>50% PB); mean = {0:4.2f} '
                     '(+{1:4.2f} - '
                      '{2:4.2f})').format(stats_int_orig_50[1],
                                          (stats_int_orig_50[3]-
                                           stats_int_orig_50[1]),
                                          (stats_int_orig_50[1] -
                                           stats_int_orig_50[2]))
             )
    
    ax1.hist(peak_ratio_orig, bins=histbins,
              histtype='step', color=colors[1],
             linestyle='solid',
             label = ('Peak flux; mean = {0:4.2f} '
                     '(+{1:4.2f} - '
                      '{2:4.2f})').format(stats_peak_orig[1],
                                          (stats_peak_orig[3]-
                                           stats_peak_orig[1]),
                                          (stats_peak_orig[1] -
                                           stats_peak_orig[2]))
             )
    

    
    ax1.hist(peak_ratio_orig_50, bins=histbins,
             histtype='step', color=colors[1],
             linestyle='dashed',
             label = ('Peak flux (>50% PB); mean = {0:4.2f} '
                     '(+{1:4.2f} - '
                      '{2:4.2f})').format(stats_peak_orig_50[1],
                                          (stats_peak_orig_50[3]-
                                           stats_peak_orig_50[1]),
                                          (stats_peak_orig_50[1] -
                                           stats_peak_orig_50[2]))
             )

    ax1.set_xlabel('Apertif / NVSS integrated flux ratio, Original GPR')
    

    ax1.set_yscale('log')
    ax1.set_ylabel('Counts')

    ax1.set_ylim([0.5,5000])

    ax1.legend()

    #med_peak_ratio = np.nanmedian(peak_ratio)
    #mean_peak_ratio = np.nanmedian(peak_ratio)
    #nd_int_16 = np.nanpercentile(peak_ratio,16)
    #nd_int_84 = np.nanpercentile(peak_ratio,84)

    #ylim = ax1.get_ylim()
    #ax1.plot([med_peak_ratio,med_peak_ratio],ylim,
    #         label = ('Peak flux ratio; '
    #                  'Median = {0:4.2f} '
    #                  "(+{1:4.2f} - {2:4.2f})").format(med_peak_ratio,
    #                                                   nd_int_16,
    #                                                   nd_int_84) )
    #ax1.legend()


    
    figpath = os.path.join(figdir,"external_flux_comp_hist.pdf")

    plt.savefig(figpath)
    plt.close()
