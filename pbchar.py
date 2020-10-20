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

#define global directories
this_dir,this_filename = os.path.split(__file__)
pbchardir = this_dir
filedir = os.path.join(pbchardir,"files")
figdir = os.path.join(pbchardir,"figures")

#create PB object
class PB(object):
    def __init__(self,
                 matchfile,
                 startdate = None,
                 enddate = None):
        """
        Initialize PB object
        Inputs:
        - matchfile (str): Path of matchfile with all cross-match info
        - startdate (str): YYMMDD, start date to make time-limited set of plots
        - enddate (str): YYMMDD, end date to make time-limited set of plots
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
            print(len(ind),len(limit_matces),len(self.matches))
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

        #make second time if date range is set
        if self.matches_date_range is not None:
            self.position_plots()
            self.oned_plots(daterange=True)

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
                       ".png").format(self.pbname,self.beam,
                                      self.startdate,self.enddate)
        else:
            peak_ratio = self.matches['peak_flux_ap']/self.matches['int_flux_nvss']
            int_ratio = self.matches['int_flux_ap']/self.matches['int_flux_nvss']
            delta_ra = self.matches['delta_ra']
            delta_dec = self.matches['delta_dec']
            pb_level = self.matches['pb_level']
            radius = self.matches['radius']
            figname = "oned_B{1}_{0}.png".format(self.pbname,self.beam)
        #setup figure; 4 plots, peak & int
        fig, axes = plt.subplots(4,2,figsize = (10,20),
                                 sharex = 'row',
                                 sharey = 'col')

        #plot function of radius
        ax1 = plot_oned_panel(axes[0,0],radius,peak_ratio)
        ax1.set_ylabel('Apertif peak flux / NVSS integrated flux')
        ax1.set_xlabel('Radius [arcmin]')
        ax1.set_title('Apertif peak / NVSS int')

        ax2 = plot_oned_panel(axes[0,1],radius,int_ratio)
        ax2.set_ylabel('Apertif int flux / NVSS integrated flux')
        ax2.set_xlabel('Radius [arcmin]')
        ax2.set_title('Apertif int / NVSS int')

        #plot function pb level
        ax3 = plot_oned_panel(axes[1,0],pb_level,peak_ratio)
        ax3.set_ylabel('Apertif peak flux / NVSS integrated flux')
        ax3.set_xlabel('Primary beam response level')
        ax3.set_xlim(1.0,0.1)

        ax4 = plot_oned_panel(axes[1,1],pb_level,int_ratio)
        ax4.set_ylabel('Apertif int flux / NVSS integrated flux')
        ax4.set_xlabel('Primary beam response level')
        ax4.set_xlim(1.0,0.1)

        #plot function ra offset
        ax5 = plot_oned_panel(axes[2,0],delta_ra,peak_ratio)
        ax5.set_ylabel('Apertif peak flux / NVSS integrated flux')
        ax5.set_xlabel('Delta RA [arcmin]')

        ax6 = plot_oned_panel(axes[2,1],delta_ra,int_ratio)
        ax6.set_ylabel('Apertif int flux / NVSS integrated flux')
        ax6.set_xlabel('Delta RA [arcmin]')

        #plot function dec offset
        ax7 = plot_oned_panel(axes[3,0],delta_dec,peak_ratio)
        ax7.set_ylabel('Apertif peak flux / NVSS integrated flux')
        ax7.set_xlabel('Delta Dec [arcmin]')

        ax8 = plot_oned_panel(axes[3,1],delta_dec,int_ratio)
        ax8.set_ylabel('Apertif int flux / NVSS integrated flux')
        ax8.set_xlabel('Delta Dec [arcmin]')
        
        #save plot to pre-defined output, based on PB/CB
        figpath = os.path.join(figdir,figname)
        plt.savefig(figpath)
        


    def peak_ratio_radius(self,daterange=False):
        """
        1-D plot of peak ratio as a function of radius
        Also add 1-D plots for delta_ra,delta_dec (offsets)
        And also as function of pb_level
        """
        #get variables, based on date range
        if daterange:
            peak_ratio = ( self.matches_date_range['peak_flux_ap'] /
                           self.matches_date_range['int_flux_nvss'] )
            delta_ra = self.matches_date_range['delta_ra']
            delta_dec = self.matches_date_range['delta_dec']
            pb_level = self.matches_date_range['pb_level']
            radius = self.matches_date_range['radius']
            figname = ("{0}_B{1}_{2}_{3}_"
                       "peak_ratio_radius.png").format(self.pbname,self.beam,
                                                         self.startdate,self.enddate)
        else:
            peak_ratio = self.matches['peak_flux_ap']/self.matches['int_flux_nvss']
            delta_ra = self.matches['delta_ra']
            delta_dec = self.matches['delta_dec']
            pb_level = self.matches['pb_level']
            radius = self.matches['radius']
            figname = "{0}_B{1}_peak_ratio_radius.png".format(self.pbname,self.beam)

        #setup figure, 4 subplots
        fig, axes = plt.subplots(2,2,figsize = (12,12),
                                 sharex = 'row',
                                 sharey = 'col')
        #plot function of radius
        axes[0,0].scatter(radius,peak_ratio,marker='.')
        axes[0,0].set_ylabel('Apertif peak flux / NVSS integrated flux')
        axes[0,0].set_xlabel('Radius [arcmin]')
        #add running mean
        #sort peak_ratio by radius
        rad, peak_mean = running_mean(radius,peak_ratio,20)
        axes[0,0].plot(rad, peak_mean,
                       label='Running mean over 20 sources',
                       color='red')
        #add relevant lines
        axes[0,0].plot([0,np.max(radius)],[1,1],color='black')
        axes[0,0].plot([0,np.max(radius)],[1.2,1.2],color='black',linestyle='--')
        axes[0,0].plot([0,np.max(radius)],[0.8,0.8],color='black',linestyle='--')
        #plot funciton of PB level
        axes[0,1].scatter(pb_level,peak_ratio,marker='.')
        axes[0,1].set_xlabel('Primary beam response level')
        axes[0,1].plot([1,0.1],[1,1],color='black')
        axes[0,1].plot([1,0.1],[1.2,1.2],color='black',linestyle='--')
        axes[0,1].plot([1,0.1],[0.8,0.8],color='black',linestyle='--')
        axes[0,1].set_xlim(1,0.1)
        lev, peak_mean = running_mean(pb_level,peak_ratio,20)
        axes[0,1].plot(lev, peak_mean,
                       label='Running mean over 20 sources',
                       color='red')
        #plot function of RA offset
        axes[1,0].scatter(delta_ra,peak_ratio,marker='.')
        axes[1,0].set_ylabel('Apertif peak flux / NVSS integrated flux')
        axes[1,0].set_xlabel('Delta R.A. [arcmin]')
        axes[1,0].plot([-30,30],[1,1],color='black')
        axes[1,0].plot([-30,30],[1.2,1.2],color='black',linestyle='--')
        axes[1,0].plot([-30,30],[0.8,0.8],color='black',linestyle='--')
        ra, peak_mean = running_mean(delta_ra,peak_ratio,20)
        axes[1,0].plot(ra, peak_mean,
                       label='Running mean over 20 sources',
                       color='red')
        #plot function of dec offset
        axes[1,1].scatter(delta_dec,peak_ratio,marker='.')
        #axes[1,1].set_ylabel('Apertif peak flux / NVSS integrated flux')
        axes[1,1].set_xlabel('Delta Decl. [arcmin]')
        axes[1,1].plot([-30,30],[1,1],color='black')
        axes[1,1].plot([-30,30],[1.2,1.2],color='black',linestyle='--')
        axes[1,1].plot([-30,30],[0.8,0.8],color='black',linestyle='--')
        dec, peak_mean = running_mean(delta_dec,peak_ratio,20)
        axes[1,1].plot(dec, peak_mean,
                       label='Running mean over 20 sources',
                       color='red')
        #save plot to pre-defined output, based on PB/CB
        figpath = os.path.join(figdir,figname)
        plt.savefig(figpath)


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

def running_mean(x,y,N):
    #cumsum = np.cumsum(np.insert(x,0,0))
    #return (cumsum[N:] - cumsum[:-N]) / float(N)
    x_inds = x.argsort()
    y_sorted = y[x_inds]
    x_sorted = x[x_inds]
    y_running_mean = np.convolve(y_sorted, np.ones((N,))/N, mode='same')
    return x_sorted,y_running_mean

def plot_oned_panel(ax,x,y):
    """
    Helper function that takes x, y arrays for a
    fig/ax pair and does the plotting
    Since I'm doing a lot of repetitive plotting
    """
    xs, ys_mean = running_mean(x,y,20)
    ax.scatter(x,y,marker='.')
    ax.plot(xs,ys_mean,label='Running mean over 20 sources',
            color='red')
    ax.plot([np.min(x),np.max(x)],[1,1],color='black')
    ax.plot([np.min(x),np.max(x)],[1.2,1.2],color='black',linestyle='--')
    ax.plot([np.min(x),np.max(x)],[0.8,0.8],color='black',linestyle='--')

    return ax
