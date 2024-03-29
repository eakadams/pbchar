#PBchar for given CB

from __future__ import print_function

__author__ = "E.A.K. Adams"

"""
Code for doing PB characterization/analysis
Create a PB object that is initiatlized with
combined cross-match file
Then run it to create different plots and such

Improvements / to-do list:
- Add optional date range to make separate set of plots
"""

from astropy.io import ascii
import numpy as np
import os
import matplotlib.pyplot as plt
from scipy.stats import norm
import matplotlib.mlab as mlab
from astropy.coordinates import SkyCoord
from astropy import units as u


#define global directories
this_dir,this_filename = os.path.split(__file__)
pbchardir = this_dir
filedir = os.path.join(pbchardir,"files")
figdir = os.path.join(pbchardir,"figures")

#get mpl colors
prop_cycle = plt.rcParams['axes.prop_cycle']
mpcolors = prop_cycle.by_key()['color']

#create PB object
class PB(object):
    def __init__(self,
                 matchfile,
                 startdate = None,
                 enddate = None,
                 N=20,
                 SN = None,
                 raterr = None,
                 size = None,
                 mode = "int"):
        """
        Initialize PB object
        Inputs:
        - matchfile (str): Path of matchfile with all cross-match info
        - startdate (str): YYMMDD, start date to make time-limited set of plots
        - enddate (str): YYMMDD, end date to make time-limited set of plots
        - N (int): number of sources to average over for binnning / running mean 
        - SN (float): optional filter on SN of measured APertif sources
        - raterr (float): optional filter on error in flux ratio
        - size (float): optional filter on NVSS size
        """
        #read in match file
        self.matches = ascii.read(matchfile,format='csv')

        #get names
        basename = os.path.basename(matchfile)
        names = basename.split("_")
        self.pbname = names[0]
        names2 = names[1].split(".")
        self.beam = names2[0]
        #print(self.pbname,self.beam)

        #add start and end date as attributes
        self.startdate = startdate
        self.enddate = enddate

        #add number of sources to use for running mean / binning
        self.N = N
        #add filters/other options
        self.SN = SN
        self.raterr = raterr
        self.size = size
        self.mode = mode

        #get a limited set of matches if date ranges are set
        #initialize attribute
        self.matches_date_range = None
        if (self.startdate is not None) or (self.enddate is not None):
            self.get_time_range_matches()

    def get_time_range_matches(self):
        """
        Get a time-range limited table of matches
        """
        #start with start date
        if self.startdate is not None:
            startid_str = self.startdate + '000'
            startid = np.int(startid_str)
            ind = np.where(self.matches['ObsID'] > startid)[0]
            print(len(ind),len(self.matches))
            limit_matches = self.matches[ind]
        else:
            limit_matches = self.matches

        #then do end date, working on limit_matches
        if self.enddate is not None:
            endid_str = self.enddate+'999'
            endid = np.int(endid_str)
            ind = np.where(limit_matches['ObsID'] < endid)[0]
            print(len(ind),len(limit_matches),len(self.matches))
            self.matches_date_range = limit_matches[ind]
        else:
            self.matches_date_range = limit_matches

    def go(self):
        """
        Run all available analysis plots:
        - position (peak & int ratios)
        - oned (for both peak & int)

        Optionally, same plots also date-range-limited:
        """
        #make plots
        self.position_plots()
        self.oned_plots()
        self.pb_level_plots()
        self.plot_hist_peak()
        
        #make second time if date range is set
        if self.matches_date_range is not None:
            self.position_plots(daterange=True)
            self.oned_plots(daterange=True)

    def internal_comp(self):
        """
        Does internal flux comparison using MDS fields

        Parameters
        ----------
        self : PB object

        Returns
        -------
        ????
        """
        #get NVSS name
        #to find unique NVSS sources
        #can't run unique on NVSS coords directly
        NVSS_name = np.empty(len(self.matches),dtype=object)
        for i in range(len(self.matches)):
            NVSS_name[i] = 'NVSS{0}+{1}'.format(self.matches['ra_nvss'][i],
                                                self.matches['dec_nvss'][i])
        #get unique NVSS sources / number of occurences
        unique_NVSS, n_counts = np.unique(NVSS_name, return_counts = True)

        #find (indices of) repeated sources
        ind_repeat = np.where(n_counts > 1)[0]
        repeated_NVSS = unique_NVSS[ind_repeat]
        
        #(iterate through repeated sources to) get "global" comparison values
        n_visits = n_counts[ind_repeat]
        median_apertif_int_flux = np.empty(len(ind_repeat))
        mean_apertif_int_flux = np.empty(len(ind_repeat))
        rms_apertif_int_flux = np.empty(len(ind_repeat))
        median_apertif_peak_flux = np.empty(len(ind_repeat))
        mean_apertif_peak_flux = np.empty(len(ind_repeat))
        rms_apertif_peak_flux = np.empty(len(ind_repeat))
        for i,source in enumerate(repeated_NVSS):
            ind_aper = np.where(NVSS_name == source)[0]
            median_apertif_int_flux[i] = np.median(
                self.matches['int_flux_ap'][ind_aper])
            mean_apertif_int_flux[i] = np.mean(
                self.matches['int_flux_ap'][ind_aper])
            rms_apertif_int_flux[i] = np.sqrt(np.mean(
                self.matches['int_flux_ap'][ind_aper]**2))
            median_apertif_peak_flux[i] = np.median(
                self.matches['peak_flux_ap'][ind_aper])
            mean_apertif_peak_flux[i] = np.mean(
                self.matches['peak_flux_ap'][ind_aper])
            rms_apertif_peak_flux[i] = np.sqrt(np.mean(
                self.matches['peak_flux_ap'][ind_aper]**2))

        #now get table of sources with individual apertif measurements
        #create a new table with added columns
        #fill columns nans and later remove nan elements

        comp_table = self.matches
        comp_table['NVSS_name'] = NVSS_name
        #add columns
        comp_table['med_ap_int_flux'] = np.full(len(comp_table),np.nan)
        comp_table['mean_ap_int_flux'] = np.full(len(comp_table),np.nan)
        comp_table['rms_ap_int_flux'] = np.full(len(comp_table),np.nan)
        comp_table['med_ap_peak_flux'] = np.full(len(comp_table),np.nan)
        comp_table['mean_ap_peak_flux'] = np.full(len(comp_table),np.nan)
        comp_table['rms_ap_peak_flux'] = np.full(len(comp_table),np.nan)
        comp_table['n_visits'] = np.full(len(comp_table),1)
        
        #iterate through and fill columns where there are multiple obs
        for i,entry in enumerate(NVSS_name):
            ind_nvss = np.where(entry == repeated_NVSS)[0]
            if len(ind_nvss) > 0:
                #entry is in repeated NVSS list
                comp_table['med_ap_int_flux'][i] = median_apertif_int_flux[ind_nvss]
                comp_table['mean_ap_int_flux'][i] = mean_apertif_int_flux[ind_nvss]
                comp_table['rms_ap_int_flux'][i] = rms_apertif_int_flux[ind_nvss]
                comp_table['med_ap_peak_flux'][i] = median_apertif_peak_flux[ind_nvss]
                comp_table['mean_ap_peak_flux'][i] = mean_apertif_peak_flux[ind_nvss]
                comp_table['rms_ap_peak_flux'][i] = rms_apertif_peak_flux[ind_nvss]
                #also get number of visits
                comp_table['n_visits'][i] = n_visits[ind_nvss]
                
        #only keep table entries w/out NaNs
        ind_repeats = np.where(comp_table['n_visits'] > 1)[0]
        comp_table = comp_table[ind_repeats]
            
        #return unique_NVSS, n_counts
        return comp_table

        
        

    def get_astrometry(self):
        """Gets astrometry for a PB object

        Parameters
        ----------
        self : PB object
        
        Returns
        -------
        offset_ra : Numpy array of Angles (astropy.coordinates.angles.Angle)
            The offset in RA to go from Apertif to NVSS coord for all matches
        offset_dec : Numpy array of Angle
            The offset in Dec to go from Apertif to NVSS coord for all matches
        median_ra_offset : float
            The median RA offset in arcsec over all observations included in PB object
        median_dec_offset : float
            The median Dec offset in arcsec over all observations in PB object
        rms_ra_offset : float
            RMS of RA offset
        rms_dec_offset : float
            RMS of Dec offset
        obsid : numpy array
            All Observation IDs in PB object
        med_ra_offset_obsid : Numpy array of Angles
            Median RA offset for each ObsID in obsid 
        med_dec_offset_obsid : Numpy array of Angles
            Median Dec offset for each ObsID in obsid 
        """

        #Load RA,Dec into skycoord objects
        coords_nvss = SkyCoord(self.matches['ra_nvss'],self.matches['dec_nvss'],
                               unit='deg')
        coords_apertif  = SkyCoord(self.matches['ra_apertif'],
                                   self.matches['dec_apertif'], unit='deg')

        #get spherical coord offsets to go from Apertif to NVSS
        offset_ra, offset_dec = coords_apertif.spherical_offsets_to(coords_nvss)

        #calculate median ra, dec offset
        #put in arcsec
        median_ra_offset = np.median(offset_ra).arcsec
        median_dec_offset = np.median(offset_dec).arcsec

        #also get std and report that
        rms_ra_offset = np.sqrt(np.mean(offset_ra.arcsec**2))
        rms_dec_offset = np.sqrt(np.mean(offset_dec.arcsec**2))

        #also look at total offset
        separation = coords_apertif.separation(coords_nvss)
        median_sep = np.median(separation.arcsec)
        rms_sep = np.sqrt(np.mean(separation.arcsec**2))
        #rms_offset = np.sqrt(np.mean(offset
        
        print(("To go from Apertif to NVSS, "
               "the median offset in RA is {0:4.2f} "
               "({2:4.2f}) arcsec. "
               "The median offset in Dec is {1:4.2f} "
               "({3:4.2f}) arcsec ").format(median_ra_offset,
                                            median_dec_offset,
                                            rms_ra_offset,
                                            rms_dec_offset))

        print(("The median total offset is {0:4.2f}arcesc "
               "and the rms of the total offset is "
               "{1:4.2f} arcsec").format(median_sep,rms_sep))

        #also do things based on ObsID
        obsid = np.unique(self.matches['ObsID'])
        med_ra_offset_obsid = np.empty(len(obsid),dtype=object)
        med_dec_offset_obsid = np.empty(len(obsid),dtype=object)
        for i,obs in enumerate(obsid):
            source_ind = np.where(self.matches['ObsID'] == obs)[0]
            offset_ra_obsid = offset_ra[source_ind]
            offset_dec_obsid = offset_dec[source_ind]
            med_ra_offset_obsid[i] = np.median(offset_ra_obsid)
            med_dec_offset_obsid[i] = np.median(offset_dec_obsid)
        

        return (offset_ra, offset_dec, median_ra_offset, median_dec_offset,
                rms_ra_offset, rms_dec_offset,
                obsid, med_ra_offset_obsid, med_dec_offset_obsid,
                median_sep,rms_sep)
    

    def plot_astrometry(self):
        """Makes plot of astrometry

        Parameters
        ----------
        self : PB object

        Returns
        -------
        fig2d : plt figure instance of 2D astrometry plot
        ax2d : ax instance of 2D astrometry plot
        fig_obsid : plt figure instance of obsid behavior
        ax_obsid : ax instance of obsid behavior
        """
        pass

    def oned_plots(self,daterange=False):
        """
        Make set of one-d plots; ratio as function of:
        -radius
        -pb level
        -delta ra
        -delta dec
        Inputs:
        - daterange (boolean): whether to use limited date range or not
        """
        #get variables, based on date range
        if daterange:
            peak_ratio = ( self.matches_date_range['peak_flux_ap'] /
                           self.matches_date_range['int_flux_nvss'] )
            int_ratio = ( self.matches_date_range['int_flux_ap'] /
                           self.matches_date_range['int_flux_nvss'] )
            delta_ra = self.matches_date_range['delta_ra']
            delta_dec = self.matches_date_range['delta_dec']
            pb_level = self.matches_date_range['pb_level']
            radius = self.matches_date_range['radius']
            figname = ("oned_B{1}_{0}_{2}_{3}"
                       ".pdf").format(self.pbname,self.beam,
                                      self.startdate,self.enddate)
            #get error bars for flux, if they exist
            #assume one exists if they all do (how I wrote code)
            if 'int_flux_ap_err' in self.matches.colnames:
                peak_ratio_err = ( peak_ratio *
                                   np.sqrt(
                                       (self.matches_date_range['peak_flux_ap_err']/
                                        self.matches_date_range['peak_flux_ap'])**2 +
                                       (self.matches_date_range['int_flux_nvss_err']/
                                        self.matches_date_range['int_flux_nvss'])**2 ) )
                int_ratio_err = ( peak_ratio *
                                   np.sqrt(
                                       (self.matches_date_range['int_flux_ap_err']/
                                        self.matches_date_range['int_flux_ap'])**2 +
                                       (self.matches_date_range['int_flux_nvss_err']/
                                        self.matches_date_range['int_flux_nvss'])**2 ) )
            else:
                peak_ratio_err = None
                int_ratio_err = None
                
        else:
            peak_ratio = self.matches['peak_flux_ap']/self.matches['int_flux_nvss']
            int_ratio = self.matches['int_flux_ap']/self.matches['int_flux_nvss']
            delta_ra = self.matches['delta_ra']
            delta_dec = self.matches['delta_dec']
            pb_level = self.matches['pb_level']
            radius = self.matches['radius']
            figname = "oned_B{1}_{0}.pdf".format(self.pbname,self.beam)
            #get error bars for flux, if they exist
            #assume one exists if they all do (how I wrote code)
            if 'int_flux_ap_err' in self.matches.colnames:
                peak_ratio_err = ( peak_ratio *
                                   np.sqrt( (self.matches['peak_flux_ap_err']/
                                             self.matches['peak_flux_ap'])**2 +
                                            (self.matches['int_flux_nvss_err']/
                                             self.matches['int_flux_nvss'])**2 ) )
                int_ratio_err = ( peak_ratio *
                                   np.sqrt( (self.matches['int_flux_ap_err']/
                                             self.matches['int_flux_ap'])**2 +
                                            (self.matches['int_flux_nvss_err']/
                                             self.matches['int_flux_nvss'])**2 ) )
            else:
                peak_ratio_err = None
                int_ratio_err = None
                
        #setup figure; 4 plots, peak & int
        fig, axes = plt.subplots(4,2,figsize = (10,20),
                                 sharex = 'row',
                                 sharey = 'col')

        #plot function of radius
        ax1 = self.plot_oned_panel(axes[0,0],radius,peak_ratio,yerr=peak_ratio_err)
        ax1.set_ylabel('Apertif peak flux / NVSS integrated flux')
        ax1.set_xlabel('Radius [arcmin]')
        ax1.set_title('Apertif peak / NVSS int')

        ax2 = self.plot_oned_panel(axes[0,1],radius,int_ratio, yerr=int_ratio_err)
        ax2.set_ylabel('Apertif int flux / NVSS integrated flux')
        ax2.set_xlabel('Radius [arcmin]')
        ax2.set_title('Apertif int / NVSS int')

        #plot function pb level
        ax3 = self.plot_oned_panel(axes[1,0],pb_level,peak_ratio,
                                   reverse=True, yerr=peak_ratio_err)
        ax3.set_ylabel('Apertif peak flux / NVSS integrated flux')
        ax3.set_xlabel('Primary beam response level')
        #ax3.set_xlim(1.0,0.1)

        ax4 = self.plot_oned_panel(axes[1,1],pb_level,int_ratio,
                                   reverse=True, yerr=int_ratio_err)
        ax4.set_ylabel('Apertif int flux / NVSS integrated flux')
        ax4.set_xlabel('Primary beam response level')
        #ax4.set_xlim(1.0,0.1)

        #plot function ra offset
        ax5 = self.plot_oned_panel(axes[2,0],delta_ra,peak_ratio,
                                   yerr = peak_ratio_err)
        ax5.set_ylabel('Apertif peak flux / NVSS integrated flux')
        ax5.set_xlabel('Delta RA [arcmin]')

        ax6 = self.plot_oned_panel(axes[2,1],delta_ra,int_ratio,
                                   yerr = int_ratio_err)
        ax6.set_ylabel('Apertif int flux / NVSS integrated flux')
        ax6.set_xlabel('Delta RA [arcmin]')

        #plot function dec offset
        ax7 = self.plot_oned_panel(axes[3,0],delta_dec,peak_ratio,
                                   yerr = peak_ratio_err)
        ax7.set_ylabel('Apertif peak flux / NVSS integrated flux')
        ax7.set_xlabel('Delta Dec [arcmin]')

        ax8 = self.plot_oned_panel(axes[3,1],delta_dec,int_ratio,
                                   yerr = int_ratio_err)
        ax8.set_ylabel('Apertif int flux / NVSS integrated flux')
        ax8.set_xlabel('Delta Dec [arcmin]')
        
        #save plot to pre-defined output, based on PB/CB
        figpath = os.path.join(figdir,figname)
        plt.savefig(figpath)

        #close figure to be safe
        plt.close('all')
        
    def plot_oned_panel(self,ax,x,y,yerr=None,
                        reverse=False,xbin=None):
        """
        Helper function that takes x, y arrays for a
        fig/ax pair and does the plotting
        Since I'm doing a lot of repetitive plotting
        Optionally reverse axes; used for PB level
        Optionally provide yerrorbars; then plot them
        """
        xs, ys_mean = running_mean(x,y,self.N)
        #lowest zorder to force points to background
        ax.errorbar(x,y,yerr=yerr,c=mpcolors[0],fmt='.',
                    linewidth=0.3,alpha=0.5, zorder=1)

        ax.plot([np.min(x),np.max(x)],[1,1],color='black')
        ax.plot([np.min(x),np.max(x)],[1.2,1.2],color='black',linestyle='--')
        ax.plot([np.min(x),np.max(x)],[0.8,0.8],color='black',linestyle='--')
        #add getting binned points and scatter
        #have on top w/ highest zorder
        xb,xbr,yb,ybsc = bin_scatter(x,y,self.N,reverse=reverse,xbin=xbin)
        ax.errorbar(xb,yb,xerr=xbr,yerr=ybsc,fmt='o',
                    color=mpcolors[1],linewidth=3,zorder =5)

        if reverse:
            #only do this for pblevel, so can set explicitly
            ax.set_xlim(1.0,0.1)

        return ax


    def position_plots(self,daterange=False):
        """
        2-D plots for int and peak flux ratios
        """
        #get variables, based on date range
        if daterange:
            peak_ratio = ( self.matches_date_range['peak_flux_ap'] /
                           self.matches_date_range['int_flux_nvss'] )
            int_ratio = ( self.matches_date_range['int_flux_ap'] /
                           self.matches_date_range['int_flux_nvss'] )
            delta_ra = self.matches_date_range['delta_ra']
            delta_dec = self.matches_date_range['delta_dec']
            figname = ("position_B{1}_{0}_{2}_{3}"
                       ".png").format(self.pbname,self.beam,
                                      self.startdate,self.enddate)
        else:
            peak_ratio = self.matches['peak_flux_ap']/self.matches['int_flux_nvss']
            int_ratio = self.matches['int_flux_ap']/self.matches['int_flux_nvss']
            delta_ra = self.matches['delta_ra']
            delta_dec = self.matches['delta_dec']
            figname = "position_B{1}_{0}.png".format(self.pbname,self.beam)

        #setup figure
        fig, (ax1,ax2) = plt.subplots(1,2,figsize=(12,6))

        sc_peak = ax1.scatter(delta_ra,delta_dec,
                              c=peak_ratio,s=20,
                              vmin = 0.75, vmax=1.5)
        fig.colorbar(sc_peak,ax=ax1)
        ax1.set_title('Apertif peak flux / NVSS integrated flux')
        ax1.set_xlabel('Delta RA [arcmin]')
        ax1.set_ylabel('Delta Dec [arcmin]')

        sc_int = ax2.scatter(delta_ra,delta_dec,
                              c=int_ratio,s=20,
                              vmin = 0.75, vmax=1.5)
        fig.colorbar(sc_int,ax=ax2)
        ax2.set_title('Apertif integrated flux / NVSS integrated flux')
        ax2.set_xlabel('Delta RA [arcmin]')
        ax2.set_ylabel('Delta Dec [arcmin]')

        figpath = os.path.join(figdir,figname)
        plt.savefig(figpath)

        plt.close('all')

    def pb_level_plots(self,daterange=False,xbins = np.arange(0.15,1.1,0.1)):
        """
        Focus on plots based on primary beam levels
        Best way to characterize things to give information to people
        """
        #get variables, based on date range
        if daterange:
            peak_ratio = ( self.matches_date_range['peak_flux_ap'] /
                           self.matches_date_range['int_flux_nvss'] )
            int_ratio = ( self.matches_date_range['int_flux_ap'] /
                           self.matches_date_range['int_flux_nvss'] )
            pb_level = self.matches_date_range['pb_level']
            figname = ("pblevel_B{1}_{0}_{2}_{3}"
                       ".pdf").format(self.pbname,self.beam,
                                      self.startdate,self.enddate)
            if 'int_flux_ap_err' in self.matches.colnames:
                peak_ratio_err = ( peak_ratio *
                                   np.sqrt( (self.matches_date_range['peak_flux_ap_err']/
                                             self.matches_date_range['peak_flux_ap'])**2 +
                                            (self.matches_date_range['int_flux_nvss_err']/
                                             self.matches_date_range['int_flux_nvss'])**2 ) )
                int_ratio_err = ( peak_ratio *
                                   np.sqrt( (self.matches_date_range['int_flux_ap_err']/
                                             self.matches_date_range['int_flux_ap'])**2 +
                                            (self.matches_date_range['int_flux_nvss_err']/
                                             self.matches_date_range['int_flux_nvss'])**2 ) )
            else:
                peak_ratio_err = None
                int_ratio_err = None
        else:
            peak_ratio = self.matches['peak_flux_ap']/self.matches['int_flux_nvss']
            int_ratio = self.matches['int_flux_ap']/self.matches['int_flux_nvss']
            pb_level = self.matches['pb_level']
            figname = "pblevel_B{1}_{0}.pdf".format(self.pbname,self.beam)
            if 'int_flux_ap_err' in self.matches.colnames:
                peak_ratio_err = ( peak_ratio *
                                   np.sqrt( (self.matches['peak_flux_ap_err']/
                                             self.matches['peak_flux_ap'])**2 +
                                            (self.matches['int_flux_nvss_err']/
                                             self.matches['int_flux_nvss'])**2 ) )
                int_ratio_err = ( peak_ratio *
                                   np.sqrt( (self.matches['int_flux_ap_err']/
                                             self.matches['int_flux_ap'])**2 +
                                            (self.matches['int_flux_nvss_err']/
                                             self.matches['int_flux_nvss'])**2 ) )
            else:
                peak_ratio_err = None
                int_ratio_err = None
                
        #setup figure; N source binning, plus equal pblevel binning, peak & int
        fig, axes = plt.subplots(2,2,figsize = (10,10),
                                 sharex = 'row',
                                 sharey = 'col')

        #plot N source binning
        ax1 = self.plot_oned_panel(axes[0,0],pb_level,peak_ratio,
                                   reverse=True, yerr= peak_ratio_err)
        ax1.set_ylabel('Apertif peak flux / NVSS integrated flux')
        ax1.set_xlabel('Primary beam response level')
        ax1.set_title('Apertif peak / NVSS int')

        ax2 = self.plot_oned_panel(axes[0,1],pb_level,int_ratio,
                                   reverse=True, yerr= int_ratio_err)
        ax2.set_ylabel('Apertif total flux / NVSS integrated flux')
        ax2.set_xlabel('Primary beam response level')
        ax2.set_title('Apertif int / NVSS int')

        #plot equal pb level binning
        
        ax3 = self.plot_oned_panel(axes[1,0],pb_level,peak_ratio,
                                   xbin=xbins,reverse=True, yerr=peak_ratio_err)
        ax3.set_ylabel('Apertif peak flux / NVSS integrated flux')
        ax3.set_xlabel('Primary beam response level')

        ax4 = self.plot_oned_panel(axes[1,1],pb_level,int_ratio,
                                   xbin=xbins,reverse=True, yerr=int_ratio_err)
        ax4.set_ylabel('Apertif int flux / NVSS integrated flux')
        ax4.set_xlabel('Primary beam response level')

        figpath = os.path.join(figdir,figname)
        plt.savefig(figpath)

        plt.close('all')

    def get_ratio_vals(self, mode='int', level = None):
        """
        Get mean, median, std and mean/median error
        Do just for integrated flux ratio; can toggle for peak
        Also for level; None is everything
        """
        #get mode working in first
        #then can filter just for that mode
        if self.mode == 'peak':
            ratio = self.matches['peak_flux_ap']/self.matches['int_flux_nvss']
            SN = self.matches['peak_flux_ap']/self.matches['peak_flux_ap_err']
            error = ( ratio *
                           np.sqrt( (self.matches['peak_flux_ap_err']/
                                     self.matches['peak_flux_ap'])**2 +
                                    (self.matches['int_flux_nvss_err']/
                                     self.matches['int_flux_nvss'])**2 ) )
        else:
            ratio = self.matches['int_flux_ap']/self.matches['int_flux_nvss']
            SN = self.matches['int_flux_ap']/self.matches['int_flux_ap_err']
            error = ( ratio *
                          np.sqrt( (self.matches['int_flux_ap_err']/
                                    self.matches['int_flux_ap'])**2 +
                                   (self.matches['int_flux_nvss_err']/
                                    self.matches['int_flux_nvss'])**2 ) )
        #independent of mode                                         
        pb_level = self.matches['pb_level']
        size = self.matches['maj_nvss']
        
        #filter on size, if set
        if self.size is not None:
            ind  = np.where(size <= self.size)[0]
            ratio = ratio[ind]
            SN = SN[ind]
            error = error[ind]
            pb_level = pb_level[ind]
            size = size[ind]
        #filter on SN, if set
        if self.SN is not None:
            ind = np.where( SN >= self.SN)[0]
            ratio = ratio[ind]
            SN = SN[ind]
            error = error[ind]
            pb_level = pb_level[ind]
            size = size[ind]
        #filter on ratio error, if set
        if self.raterr is not None:
            ind = np.where(error <= self.raterr)[0]
            ratio = ratio[ind]
            SN = SN[ind]
            error = error[ind]
            pb_level = pb_level[ind]
            size = size[ind]

        #if pb level limits
        if level is not None:
            ind = np.where(pb_level >= level)[0]
            ratio =ratio[ind]
            error = error[ind]

        #get stats
        mean_ratio = np.mean(ratio)
        median_ratio = np.median(ratio)
        std_ratio = np.std(ratio)
        mean_error =np.mean(error)
        median_error = np.median(error)

        return mean_ratio,median_ratio,std_ratio,mean_error,median_error

    def plot_hist_peak(self,daterange=False):
        """
        Plot histograms of flux ratios
        If daterange is set, do that as another histogram in same plot
        Not turned on yet (lazy, not really using this feature)
        Would be better to set date range globally and just update matches
        Add more plots - peak to see about size bias
        And filtering by NVSS size
        Also limit bins to 0.5-1.5 to avoid outliers
        Need to add helper function for histogram plots to clean code up
        """
        peak_ratio = self.matches['peak_flux_ap']/self.matches['int_flux_nvss']
        int_ratio = self.matches['int_flux_ap']/self.matches['int_flux_nvss']
        pb_level = self.matches['pb_level']
        figname = "flux_ratio_hist_B{1}_{0}.pdf".format(self.pbname,self.beam)

        fig, ((ax1,ax2),(ax3,ax4)) = plt.subplots(2,2,figsize = (10,10))

        #do integrated flux ratio
        #do norm fitting
        (mu, sigma) = norm.fit(int_ratio)
        #normalize histogram for fitting
        n, bins, patches = ax1.hist(int_ratio,
                                    label='All PB levels',
                                    density=True,
                                    bins = np.arange(0.5,1.55,0.05),
                                    alpha=0.9,
                                    histtype='step')
        ax1.set_xlabel('Apertif int flux / NVSS int flux')
        # add a 'best fit' line
        y = mlab.normpdf( bins, mu, sigma)
        ax1.plot(bins, y, linestyle='--',
                 label='Gaussian mean={0:4.2f}, sigma={1:4.2f}'.format(mu,sigma))
        
        #get only where PB level is >= 0.5
        ind_main = np.where(pb_level >= 0.5)[0]
        (mu_50, sigma_50) = norm.fit(int_ratio[ind_main])
        n, bins, patches =ax1.hist(int_ratio[ind_main],
                                   label='PB level >= 0.5',
                                   density=True,
                                   bins = np.arange(0.5,1.55,0.05),
                                   alpha=0.9,
                                   histtype='step')
        # add a 'best fit' line
        y = mlab.normpdf( bins, mu_50, sigma_50)
        ax1.plot(bins, y, linestyle='--',
                 label=('PB>50%; Gaussian mean={0:4.2f}, '
                        'sigma={1:4.2f}').format(mu_50,sigma_50))

        ax1.set_ylim(0,4.5)
        ax1.set_xlim(0.2,2.5)
        
        ax1.legend()

        #do peak flux ratio
        #do norm fitting
        (mu, sigma) = norm.fit(peak_ratio)
        #normalize histogram for fitting
        n, bins, patches = ax2.hist(peak_ratio,
                                    label='All PB levels',
                                    density=True,
                                    bins = np.arange(0.5,1.55,0.05),
                                    alpha=0.9,
                                    histtype='step')
        ax2.set_xlabel('Apertif peak flux / NVSS int flux')
        # add a 'best fit' line
        y = mlab.normpdf( bins, mu, sigma)
        ax2.plot(bins, y, linestyle='--',
                 label='Gaussian mean={0:4.2f}, sigma={1:4.2f}'.format(mu,sigma))
        
        #get only where PB level is >= 0.5
        ind_main = np.where(pb_level >= 0.5)[0]
        (mu_50, sigma_50) = norm.fit(peak_ratio[ind_main])
        n, bins, patches =ax2.hist(peak_ratio[ind_main],
                                   label='PB level >= 0.5',
                                   density=True,
                                   bins = np.arange(0.5,1.55,0.05),
                                   alpha=0.9,
                                   histtype='step')
        # add a 'best fit' line
        y = mlab.normpdf( bins, mu_50, sigma_50)
        ax2.plot(bins, y, linestyle='--',
                 label=('PB>50%; Gaussian mean={0:4.2f}, '
                        'sigma={1:4.2f}').format(mu_50,sigma_50))

        ax2.set_ylim(0,4.5)
        ax2.set_xlim(0.2,2.5)
        
        ax2.legend()


        #do integrated flux ratio w/ limited size
        ind40 = np.where(self.matches['maj_nvss'] <= 40.)[0]
        int_ratio_40 = int_ratio[ind40]
        pb_level_40 = pb_level[ind40]
        (mu, sigma) = norm.fit(int_ratio_40)
        n, bins, patches = ax3.hist(int_ratio_40,
                                    label='All PB levels',
                                    density=True,
                                    bins = np.arange(0.5,1.55,0.05),
                                    alpha=0.9,
                                    histtype='step')
        ax3.set_xlabel('Apertif int flux / NVSS int flux')
        ax3.set_title('NVSS major axis <= 40"')
        # add a 'best fit' line
        y = mlab.normpdf( bins, mu, sigma)
        ax3.plot(bins, y, linestyle='--',
                 label='Gaussian mean={0:4.2f}, sigma={1:4.2f}'.format(mu,sigma))
        
        #get only where PB level is >= 0.5
        ind_main = np.where(pb_level_40 >= 0.5)[0]
        (mu_50, sigma_50) = norm.fit(int_ratio_40[ind_main])
        n, bins, patches =ax3.hist(int_ratio_40[ind_main],
                                   label='PB level >= 0.5',
                                   density=True,
                                   bins = np.arange(0.5,1.55,0.05),
                                   alpha=0.9,
                                   histtype='step')
        # add a 'best fit' line
        y = mlab.normpdf( bins, mu_50, sigma_50)
        ax3.plot(bins, y, linestyle='--',
                 label=('PB>50%; Gaussian mean={0:4.2f}, '
                        'sigma={1:4.2f}').format(mu_50,sigma_50))
        ax3.set_ylim(0,4.5)
        ax3.set_xlim(0.2,2.5)
        ax3.legend()

        #do peak flux ratio w/ limited size
        ind40 = np.where(self.matches['maj_nvss'] <= 40.)[0]
        peak_ratio_40 = peak_ratio[ind40]
        pb_level_40 = pb_level[ind40]
        (mu, sigma) = norm.fit(peak_ratio_40)
        n, bins, patches = ax4.hist(peak_ratio_40,
                                    label='All PB levels',
                                    density=True,
                                    bins = np.arange(0.5,1.55,0.05),
                                    alpha=0.9,
                                    histtype='step')
        ax4.set_xlabel('Apertif peak flux / NVSS int flux')
        ax4.set_title('NVSS major axis <= 40"')
        # add a 'best fit' line
        y = mlab.normpdf( bins, mu, sigma)
        ax4.plot(bins, y, linestyle='--',
                 label='Gaussian mean={0:4.2f}, sigma={1:4.2f}'.format(mu,sigma))
        
        #get only where PB level is >= 0.5
        ind_main = np.where(pb_level_40 >= 0.5)[0]
        (mu_50, sigma_50) = norm.fit(peak_ratio_40[ind_main])
        n, bins, patches =ax4.hist(peak_ratio_40[ind_main],
                                   label='PB level >= 0.5',
                                   density=True,
                                   bins = np.arange(0.5,1.55,0.05),
                                   alpha=0.9,
                                   histtype='step')
        # add a 'best fit' line
        y = mlab.normpdf( bins, mu_50, sigma_50)
        ax4.plot(bins, y, linestyle='--',
                 label=('PB>50%; Gaussian mean={0:4.2f}, '
                        'sigma={1:4.2f}').format(mu_50,sigma_50))
        ax4.set_ylim(0,4.5)
        ax4.set_xlim(0.2,2.5)
        ax4.legend()

        

        figpath = os.path.join(figdir,figname)
        plt.savefig(figpath)

        plt.close('all')
    
    def get_scatter_pblev(self,lev,daterange=False,
                          xbins=np.arange(0.15,1.1,0.1),
                          mode='int'):
        """
        Get scatter as a function of lev
        Want to get scatter both for equally spaced bins 
        and bins by N sources
        Can specify bins
        Also whether do for more 'int' or 'peak' flux
        Default int
        """
        #get variables, based on date range
        if daterange:
            peak_ratio = ( self.matches_date_range['peak_flux_ap'] /
                           self.matches_date_range['int_flux_nvss'] )
            int_ratio = ( self.matches_date_range['int_flux_ap'] /
                           self.matches_date_range['int_flux_nvss'] )
            pb_level = self.matches_date_range['pb_level']
        else:
            peak_ratio = self.matches['peak_flux_ap']/self.matches['int_flux_nvss']
            int_ratio = self.matches['int_flux_ap']/self.matches['int_flux_nvss']
            pb_level = self.matches['pb_level']

        if mode == 'peak':
            ratio = peak_ratio
        else:
            ratio = int_ratio

        #get various scatter values for given lev
        #xbins = np.arange(0.1,1.1,0.1)
        sc_N_lev, sc_N_mean = scatter_val(pb_level,ratio,self.N,lev,
                                          xbin=None,reverse=True)
        sc_bin_lev, sc_bin_mean = scatter_val(pb_level,ratio,self.N,lev,
                                              xbin=xbins)

        return sc_N_lev, sc_bin_lev, sc_N_mean, sc_bin_mean

    def get_mean_ratio(self,mode='int',daterange=False):
        """
        Get mean value of ratio, for integrated or peak flux
        """
         #get variables, based on date range
        if daterange:
            peak_ratio = ( self.matches_date_range['peak_flux_ap'] /
                           self.matches_date_range['int_flux_nvss'] )
            int_ratio = ( self.matches_date_range['int_flux_ap'] /
                           self.matches_date_range['int_flux_nvss'] )
        else:
            peak_ratio = self.matches['peak_flux_ap']/self.matches['int_flux_nvss']
            int_ratio = self.matches['int_flux_ap']/self.matches['int_flux_nvss']

        if mode == 'peak':
            ratio = peak_ratio
        else:
            ratio = int_ratio

        #also try to account for standard spectral index
        #apertif central freq: 1362.5 MHz.
        #nvss cen freq:
        #can't easily find, skip this for now
        #Can update for it later
            
        return np.mean(ratio)

def scatter_val(x,y,N,val,xbin=None,reverse=False):
    """
    Get scatter for given val (x-array)
    Should also return mean value
    """
    xb,xbr,yb,ybsc = bin_scatter(x,y,N,xbin=xbin,reverse=reverse)
    xlow = np.around(xb-xbr,4)
    xhigh = np.around(xb+xbr,4)
    ind1 = np.where(xlow <= val)[0]
    ind2 = np.where(xhigh > val)[0]
    ind = np.intersect1d(ind1,ind2)
    if len(ind) == 0:
        #this happens when bin edges are close to value
        #due to floating point / rounding errors
        #find the bin edge that is closest and use it
        test1 = np.min(ind1)
        test2 = np.max(ind2)
        lowval = xlow[test1]
        highval = xhigh[test2]
        if abs(highval - val) > abs(lowval-val):
            sc_lev = ybsc[test1]
        else:
            sc_lev = ybsc[test2]
    else:
        sc_lev = ybsc[ind][0]

    mean_sc = np.mean(ybsc)

    return sc_lev, mean_sc
    

def running_mean(x,y,N):
    #cumsum = np.cumsum(np.insert(x,0,0))
    #return (cumsum[N:] - cumsum[:-N]) / float(N)
    x_inds = x.argsort()
    y_sorted = y[x_inds]
    x_sorted = x[x_inds]
    y_running_mean = np.convolve(y_sorted, np.ones((N,))/N, mode='same')
    return x_sorted,y_running_mean

def bin_scatter(x,y,N,reverse=False,xbin=None):
    """
    Take matched x,y arrays and bin by N sources
    Return xb,yb arrays of binned values
    Plus xbr, range that binned array covers
    And ybsc, scatter in the y-direction w/in bin
    If set, reverse axes in sorting
    xbin can be used to set xbins for other type of binning
    """
    #first sort by x-array
    if reverse:
        x_inds = x.argsort()[::-1]
    else:
        x_inds = x.argsort()
    y_sorted = y[x_inds]
    x_sorted = x[x_inds]
    #either bin by Nsources or equal xbins
    if xbin is None:
        #bin by  Nsources together
        #do the binning by reshaping array
        #drop last bin if needed
        x_bin_shape = x_sorted[:(x_sorted.size // N) * N].reshape(-1, N)
        y_bin_shape = y_sorted[:(y_sorted.size // N) * N].reshape(-1, N)
        yb = y_bin_shape.mean(axis=1)
        ybsc = y_bin_shape.std(axis=1)
        xbr = (x_bin_shape.max(axis=1) - x_bin_shape.min(axis=1))/2
        xb = x_bin_shape.min(axis=1) + xbr
    else:
        #do histogram binning and find scatter
        n, bins = np.histogram(x, bins=xbin)
        sy, bins = np.histogram(x, bins=xbin, weights=y)
        sy2, bins = np.histogram(x, bins=xbin, weights=y*y)
        mean = sy / n
        ybsc = np.sqrt(sy2/n - mean*mean)
        yb = mean
        #find midpoint
        xb = (bins[1:] + bins[:-1])/2
        #get range
        xbr = (bins[1:] - bins[:-1])/2
        

    return xb,xbr,yb,ybsc
