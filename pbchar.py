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

#get mpl colors
prop_cycle = plt.rcParams['axes.prop_cycle']
mpcolors = prop_cycle.by_key()['color']

#create PB object
class PB(object):
    def __init__(self,
                 matchfile,
                 startdate = None,
                 enddate = None,
                 N=20):
        """
        Initialize PB object
        Inputs:
        - matchfile (str): Path of matchfile with all cross-match info
        - startdate (str): YYMMDD, start date to make time-limited set of plots
        - enddate (str): YYMMDD, end date to make time-limited set of plots
        - N (int): number of sources to average over for binnning / running mean 
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

        #make second time if date range is set
        if self.matches_date_range is not None:
            self.position_plots(daterange=True)
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
                       ".pdf").format(self.pbname,self.beam,
                                      self.startdate,self.enddate)
        else:
            peak_ratio = self.matches['peak_flux_ap']/self.matches['int_flux_nvss']
            int_ratio = self.matches['int_flux_ap']/self.matches['int_flux_nvss']
            delta_ra = self.matches['delta_ra']
            delta_dec = self.matches['delta_dec']
            pb_level = self.matches['pb_level']
            radius = self.matches['radius']
            figname = "oned_B{1}_{0}.pdf".format(self.pbname,self.beam)
        #setup figure; 4 plots, peak & int
        fig, axes = plt.subplots(4,2,figsize = (10,20),
                                 sharex = 'row',
                                 sharey = 'col')

        #plot function of radius
        ax1 = self.plot_oned_panel(axes[0,0],radius,peak_ratio)
        ax1.set_ylabel('Apertif peak flux / NVSS integrated flux')
        ax1.set_xlabel('Radius [arcmin]')
        ax1.set_title('Apertif peak / NVSS int')

        ax2 = self.plot_oned_panel(axes[0,1],radius,int_ratio)
        ax2.set_ylabel('Apertif int flux / NVSS integrated flux')
        ax2.set_xlabel('Radius [arcmin]')
        ax2.set_title('Apertif int / NVSS int')

        #plot function pb level
        ax3 = self.plot_oned_panel(axes[1,0],pb_level,peak_ratio, reverse=True)
        ax3.set_ylabel('Apertif peak flux / NVSS integrated flux')
        ax3.set_xlabel('Primary beam response level')
        #ax3.set_xlim(1.0,0.1)

        ax4 = self.plot_oned_panel(axes[1,1],pb_level,int_ratio, reverse=True)
        ax4.set_ylabel('Apertif int flux / NVSS integrated flux')
        ax4.set_xlabel('Primary beam response level')
        #ax4.set_xlim(1.0,0.1)

        #plot function ra offset
        ax5 = self.plot_oned_panel(axes[2,0],delta_ra,peak_ratio)
        ax5.set_ylabel('Apertif peak flux / NVSS integrated flux')
        ax5.set_xlabel('Delta RA [arcmin]')

        ax6 = self.plot_oned_panel(axes[2,1],delta_ra,int_ratio)
        ax6.set_ylabel('Apertif int flux / NVSS integrated flux')
        ax6.set_xlabel('Delta RA [arcmin]')

        #plot function dec offset
        ax7 = self.plot_oned_panel(axes[3,0],delta_dec,peak_ratio)
        ax7.set_ylabel('Apertif peak flux / NVSS integrated flux')
        ax7.set_xlabel('Delta Dec [arcmin]')

        ax8 = self.plot_oned_panel(axes[3,1],delta_dec,int_ratio)
        ax8.set_ylabel('Apertif int flux / NVSS integrated flux')
        ax8.set_xlabel('Delta Dec [arcmin]')
        
        #save plot to pre-defined output, based on PB/CB
        figpath = os.path.join(figdir,figname)
        plt.savefig(figpath)

        #close figure to be safe
        plt.close('all')
        
    def plot_oned_panel(self,ax,x,y,reverse=False,xbin=None):
        """
        Helper function that takes x, y arrays for a
        fig/ax pair and does the plotting
        Since I'm doing a lot of repetitive plotting
        Optionally reverse axes; used for PB level
        """
        xs, ys_mean = running_mean(x,y,self.N)
        ax.scatter(x,y,marker='.',c=mpcolors[0], s=5)
        #ax.plot(xs,ys_mean,
        #        label='Running mean over {} sources'.format(self.N),
        #        color=mpcolors[1])
        ax.plot([np.min(x),np.max(x)],[1,1],color='black')
        ax.plot([np.min(x),np.max(x)],[1.2,1.2],color='black',linestyle='--')
        ax.plot([np.min(x),np.max(x)],[0.8,0.8],color='black',linestyle='--')
        #add getting binned points and scatter
        xb,xbr,yb,ybsc = bin_scatter(x,y,self.N,reverse=reverse,xbin=xbin)
        ax.errorbar(xb,yb,xerr=xbr,yerr=ybsc,fmt='.',
                    color=mpcolors[1])

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

    def pb_level_plots(self,daterange=False):
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
        else:
            peak_ratio = self.matches['peak_flux_ap']/self.matches['int_flux_nvss']
            int_ratio = self.matches['int_flux_ap']/self.matches['int_flux_nvss']
            pb_level = self.matches['pb_level']
            figname = "pblevel_B{1}_{0}.pdf".format(self.pbname,self.beam)
        #setup figure; N source binning, plus equal pblevel binning, peak & int
        fig, axes = plt.subplots(2,2,figsize = (10,10),
                                 sharex = 'row',
                                 sharey = 'col')

        #plot N source binning
        ax1 = self.plot_oned_panel(axes[0,0],pb_level,peak_ratio,reverse=True)
        ax1.set_ylabel('Apertif peak flux / NVSS integrated flux')
        ax1.set_xlabel('Primary beam response level')
        ax1.set_title('Apertif peak / NVSS int')

        ax2 = self.plot_oned_panel(axes[0,1],pb_level,int_ratio,reverse=True)
        ax2.set_ylabel('Apertif total flux / NVSS integrated flux')
        ax2.set_xlabel('Primary beam response level')
        ax2.set_title('Apertif int / NVSS int')

        #plot equal pb level binning
        xbins = np.arange(0.1,1.1,0.1)
        print(xbins)
        ax3 = self.plot_oned_panel(axes[1,0],pb_level,peak_ratio,
                                   xbin=xbins,reverse=True)
        ax3.set_ylabel('Apertif peak flux / NVSS integrated flux')
        ax3.set_xlabel('Primary beam response level')

        ax4 = self.plot_oned_panel(axes[1,1],pb_level,int_ratio,
                                   xbin=xbins,reverse=True)
        ax4.set_ylabel('Apertif int flux / NVSS integrated flux')
        ax4.set_xlabel('Primary beam response level')

        figpath = os.path.join(figdir,figname)
        plt.savefig(figpath)

        plt.close('all')
        

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
        print(bins)
        print(len(n),len(bins))
        #find midpoint
        xb = (bins[1:] + bins[:-1])/2
        #get range
        xbr = (bins[1:] - bins[:-1])/2
        

    return xb,xbr,yb,ybsc
