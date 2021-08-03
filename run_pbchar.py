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
             ('mean_error','f8'),('median_error','f8'),
             ('median_ra_offset','f8'),('rms_ra_offset','f8'),
             ('median_dec_offset','f8'),('rms_dec_offset','f8'),
             ('median_total_offset','f8'),('rms_total_offset','f8')]
    #two tables; one all vals, one 50%
    t = Table(data = np.zeros(N,dtype=dtype))
    t50 = Table(data = np.zeros(N,dtype=dtype))

    #get a table with values for 2-D position plot
    #is this where I want to record array of ra/dec offsets for 2d plot?
    #how do I know that positions line up with the delta_ra/dec?
    #need to confirm / test this somehow
    dtype_pos = [('beam','i4'),('delta_ra','f8'),
                 ('delta_ra_cb','f8'), ('delta_dec','f8'),
                 ('delta_dec_cb','f8'), ('int_flux_ratio','f8'),
                 ('ra_offset','f8'),
                 ('dec_offset','f8')]
    t_pos = Table(data=np.zeros(N,dtype=dtype_pos))
    
    #run in serial per beam
    #weird issues with matplotlib / backend
    #and this is fast
    #also, want to create a single table for internal comparison
    #thus appending things
    tablelist = []
    for bm in range(40):
        #get file
        matchfile = os.path.join(args.path,"{0}_{1:02d}.csv".format(
            args.matchfilebase,bm))
        #load PB object
        CB = pbchar.PB(matchfile,startdate=args.startdate,
                       enddate=args.enddate, SN=args.SN,
                       size = args.size, raterr=args.raterr,
                       mode = args.mode)

        #get internal flux comparison tables
        comp_table = CB.internal_comp()
        tablelist.append(comp_table)
        
        #get plots
        #print(args.plots,type(args.plots))
        if args.plots:
            CB.go()
        #get scatter values for table
        (mean_ratio,median_ratio,std_ratio,
         mean_error,median_error) = CB.get_ratio_vals()
        (mean_ratio_50, median_ratio_50, std_ratio_50,
         mean_error_50, median_error_50) = CB.get_ratio_vals(level=0.5)

        (offset_ra, offset_dec, median_ra_offset,
         median_dec_offset, rms_ra_offset, rms_dec_offset,
         obsid, med_ra_offset_obsid,
         med_dec_offset_obsid, median_total_offset,
         rms_total_offset) = CB.get_astrometry()

        #can't run this at 50% level; want to actually set this in object init
        #and then in running ofthis code
        #place holders to get code to work
        median_ra_offset_50 = median_ra_offset
        median_dec_offset_50 = median_dec_offset
        rms_ra_offset_50 = rms_ra_offset
        rms_dec_offset_50 = rms_dec_offset
        median_total_offset_50 = median_total_offset
        rms_total_offset_50 = rms_total_offset
        
        #set trable row values
        t[bm] = (bm,mean_ratio,median_ratio,std_ratio,
                 mean_error,median_error, median_ra_offset,
                 rms_ra_offset,
                 median_dec_offset, rms_dec_offset,
                 median_total_offset,rms_total_offset)
        t50[bm] = (bm,mean_ratio_50, median_ratio_50,
                   std_ratio_50, mean_error_50, median_error_50,
                   median_ra_offset_50, rms_ra_offset_50,
                   median_dec_offset_50, rms_dec_offset_50,
                   median_total_offset_50, rms_total_offset_50)

    #exit work over beams and collect information

    
        
    #get a limited table for docs and write that out
    t_docs = t['beam','median_flux_ratio','std_flux_ratio','median_error',
               'median_ra_offset','rms_ra_offset','median_dec_offset',
               'rms_dec_offset','median_total_offset','rms_total_offset']
    t50_docs = t50['beam','median_flux_ratio','std_flux_ratio','median_error']
    tablename = "ratio_values_doc_{}.csv".format(args.matchfilebase)
    tablepath = os.path.join(args.path,tablename)
    t_docs.write(tablepath,format='csv',overwrite=True,
                 formats = {'median_flux_ratio': '4.2f',
                            'std_flux_ratio': '4.2f',
                            'median_error':'4.2f',
                            'median_ra_offset':'5.3f',
                            'rms_ra_offset':'5.3f',
                            'median_dec_offset':'5.3f',
                            'rms_dec_offset':'5.3f'})
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

    #get median of median offset values
    print(("The median RA offset across the median of "
           "all beams is {0:5.3f} ({1:5.3f})  arcsec").format(
               np.median(t['median_ra_offset']),
               np.sqrt(np.mean(t['rms_ra_offset']**2))))
    print(("The median Dec offset across the median of "
           "all beams is {0:5.3f} ({1:5.3f})  arcsec").format(
               np.median(t['median_dec_offset']),
               np.sqrt(np.mean(t['rms_dec_offset']**2))))

    print(("The median total offset across the median of "
           "all beams is {0:5.3f} ({1:5.3f}) arcsec").format(
               np.median(t['median_total_offset']),
               np.sqrt(np.mean(t['rms_total_offset']**2))))
          

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


    #make astrometry figures
    #want to show median values per compound beam
    #just as a check for systematics
    #want a single panel figure
    #arrow at location of each beam showing median offset
    #multiply by 100x to show up on scale?
    #color arrows by median or rms of total offset?

    #start with cbpos file
    cb_pos = ascii.read(os.path.join(filedir,"cb_offsets.txt"))
    
    #create arrays for plotting arrows
    xvals = cb_pos['ra']
    yvals = cb_pos['dec']
    #arrow goes to x+dx
    #offsets are measured from APertif to NVSS
    #want arrow showing direction from NVSS to Apertif,
    #so reverse the sign
    #mutliply by 100x and put in degrees
    delta_x = xvals - t['median_ra_offset'] * 100/3600.
    delta_y = yvals - t['median_dec_offset'] * 100/3600.
    #get color array
    colors = t['median_total_offset']
    
    fig, (ax1) = plt.subplots(1,1,figsize=(6,5))

    Q = ax1.quiver(xvals, yvals, delta_x, delta_y, colors)

    ax1.set_xlim(1.8,-1.8)
    ax1.set_ylim(-1.4,1.4)

    ax1.set_xlabel('RA [deg]')
    ax1.set_ylabel('Decl. [deg]')
    #add colorbar
    cb = fig.colorbar(Q)
    cb.set_label('Median total offset [arcsec]')

    if args.output is None:
        figpath = os.path.join(figdir,
                               "astrometry_{}.pdf".format(args.matchfilebase))
    else:
        figpath = os.path.join(figdir,
                               "astrometry_{0}_{1}.pdf".format(
                                   args.matchfilebase,
                                   args.output))
    plt.savefig(figpath)
    plt.close()

    #create a vstacked comparison table
    comp_table_all_beams = vstack(tablelist)
    print(comp_table_all_beams.colnames)
    print(len(comp_table_all_beams))

    #do filtering if set in command
    if args.SN is not None:
        #have I set SN to filter on Apertif and/or NVSS?
        #let's do both
        #first nvss
        nvss_sn = comp_table_all_beams['int_flux_nvss'] / comp_table_all_beams['int_flux_nvss_err']
        ind_nvss_sn = np.where(nvss_sn >= args.SN)[0]
        comp_table_all_beams = comp_table_all_beams[ind_nvss_sn]
        #then apertif, both int and peak
        ap_int_sn = comp_table_all_beams['int_flux_ap'] / comp_table_all_beams['int_flux_ap_err']
        ind_ap_int_sn = np.where(ap_int_sn >= args.SN)[0]
        comp_table_all_beams = comp_table_all_beams[ind_ap_int_sn]
        ap_peak_sn = comp_table_all_beams['peak_flux_ap'] / comp_table_all_beams['peak_flux_ap_err']
        ind_ap_peak_sn = np.where(ap_peak_sn >= args.SN)[0]
        comp_table_all_beams = comp_table_all_beams[ind_ap_peak_sn]

    #and do size filtering. this is clearly set on NVSS size
    if args.size is not None:
        ind_size = np.where(comp_table_all_beams['maj_nvss'] <= args.size)[0]
        comp_table_all_beams = comp_table_all_beams[ind_size]

    #after all filtering check length again
    print(len(comp_table_all_beams))

    #get normalized flux differences
    normdiff_int_flux_mean = ( ( comp_table_all_beams['int_flux_ap'] -
                                 comp_table_all_beams['mean_ap_int_flux'] ) /
                               comp_table_all_beams['mean_ap_int_flux'] ) 
    normdiff_peak_flux_mean = ( ( comp_table_all_beams['peak_flux_ap'] -
                                  comp_table_all_beams['mean_ap_peak_flux'] ) /
                                comp_table_all_beams['mean_ap_peak_flux'] )
    
    #make figures for internal flux comparison
    #start with integrated; reference the Apertif median
    #color by number of measurements

    fig, ((ax1, ax2) )= plt.subplots(1,2,figsize=(11,5))

    #get colors based on nvisits
    #require more than two visits
    #see a huge dichotomy that has to do with fields with very different
    #flux scaling. Probably observations where something just went bad

    ind_3 = np.where(comp_table_all_beams['n_visits'] > 2)[0]

    """
    sc = ax1.scatter(comp_table_all_beams['med_ap_int_flux'][ind_3]*1000.,
                     (((comp_table_all_beams['int_flux_ap'][ind_3] -
                        comp_table_all_beams['med_ap_int_flux'][ind_3]) /
                       comp_table_all_beams['med_ap_int_flux'][ind_3]) ),
                     marker='.',
                     c = comp_table_all_beams['n_visits'][ind_3]
    )

    sc2 = ax2.scatter(comp_table_all_beams['med_ap_peak_flux'][ind_3]*1000.,
                      (((comp_table_all_beams['peak_flux_ap'][ind_3] -
                         comp_table_all_beams['med_ap_peak_flux'][ind_3]) /
                        comp_table_all_beams['med_ap_peak_flux'][ind_3]) ),
                      marker='.',
                      c = comp_table_all_beams['n_visits'][ind_3]
    )
    """
    sc3 = ax1.scatter(comp_table_all_beams['mean_ap_int_flux'][ind_3]*1000.,
                      normdiff_int_flux_mean[ind_3],
                      marker='.', s=5,
                      c = comp_table_all_beams['n_visits'][ind_3]
    )

    sc4 = ax2.scatter(comp_table_all_beams['mean_ap_peak_flux'][ind_3]*1000.,
                      normdiff_peak_flux_mean[ind_3],
                      marker='.', s=5,
                      c = comp_table_all_beams['n_visits'][ind_3]
    )

    cb = fig.colorbar(sc3)

    """
    ax1.set_xscale('log')
    #ax1.set_yscale('log')
    ax1.set_xlabel('median Apertif integrated flux [mJy]')
    ax1.set_ylabel('(Individual - median) / median integrated flux')

    ax2.set_xscale('log')
    ax2.set_xlabel('median Apertif peak flux [mJy]')
    ax2.set_ylabel('(Individual - median) / median peak flux')
    """
    
    ax1.set_xscale('log')
    ax1.set_xlabel('mean Apertif integrated flux [mJy]')
    ax1.set_ylabel('(Individual - mean) / mean integrated flux')

    ax2.set_xscale('log')
    ax2.set_xlabel('mean Apertif peak flux [mJy]')
    ax2.set_ylabel('(Individual - mean) / mean peak flux')
        

    if args.output is None:
        figpath = os.path.join(figdir,
                               "internal_flux_comp_{}.pdf".format(args.matchfilebase))
    else:
        figpath = os.path.join(figdir,
                               "internal_flux_comp_{0}_{1}.pdf".format(
                                   args.matchfilebase,
                                   args.output))
    plt.savefig(figpath)
    plt.close()


    #add a histogram of normalized flux accuracies

    fig, ((ax1, ax2) )= plt.subplots(1,2,figsize=(10,5))#,
                                     #sharex=True, sharey=True)

    n_bins = 50

    ax1.hist(normdiff_int_flux_mean[ind_3],
             bins = np.arange(-1.5,1.55,0.05),
             histtype='step')
    ax1.set_xlabel('(individual - mean) / mean integrated flux')
    

    ax1.set_yscale('log')
    ax1.set_ylabel('Counts')
    #ax1.set_ylim(5e-4,10)

    med_normdiff_int = np.nanmedian(normdiff_int_flux_mean[ind_3])
    #and get one sigma range
    nd_int_16 = np.nanpercentile(normdiff_int_flux_mean[ind_3],16)
    nd_int_84 = np.nanpercentile(normdiff_int_flux_mean[ind_3],84)

    ylim = ax1.get_ylim()
    ax1.plot([med_normdiff_int,med_normdiff_int],ylim,
             label = ('Normalized integrated flux; '
                      'Median = {0:4.2f} '
                      "(+{1:4.2f} - {2:4.2f})").format(med_normdiff_int,
                                                       nd_int_16,
                                                       nd_int_84) )
    ax1.legend()
    
    #(mu_int, sigma_int) = norm.fit(normdiff_int_flux_mean[ind_3])

    #print(mu_int, sigma_int)

    #y = mlab.normpdf( bins, mu_int, sigma_int)
    #l = ax1.plot(bins, y)
    
    ax2.hist( normdiff_peak_flux_mean[ind_3],
              bins = n_bins, density = True
              )

    (mu_peak, sigma_peak) = norm.fit(normdiff_peak_flux_mean[ind_3])

    print(mu_peak, sigma_peak)

    #y = mlab.normpdf( bins, mu_peak, sigma_peak)
    #l = ax2.plot(bins, y)
    ax2.set_xlabel('(individual - mean) / mean peak flux')
    

    ax2.set_yscale('log')
    ax2.set_ylim(5e-4,10)
    
    if args.output is None:
        figpath = os.path.join(figdir,
                               "internal_flux_comp_hist_{}.pdf".format(args.matchfilebase))
    else:
        figpath = os.path.join(figdir,
                               "internal_flux_comp_hist_{0}_{1}.pdf".format(
                                   args.matchfilebase,
                                   args.output))
    plt.savefig(figpath)
    plt.close()
