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

#define global directories
this_dir,this_filename = os.path.split(__file__)
pbchardir = this_dir
filedir = os.path.join(pbchardir,"files")
figdir = os.path.join(pbchardir,"figures")

#create PB object
class PB(object):
    def __init__(self,
                 matchfile,
                 startdate = None
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
        print(self.pbname,self.beam)

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
        - peak_ratio_position
        - int_ratio_position

        Optionally, same plots also date-range-limited:
        - peak_ratio_position_date_range
        - int_ratio_position
        """
        #make plots
        self.peak_ratio_position()
        self.int_ratio_position()

        #make second time if date range is set
        if self.matches_date_range is not None:
            self.peak_ratio_position(daterange=True)
            self.int_ratio_position(daterange=True)


    def peak_ratio_position(self,daterange=False):
        """
        Plot peak ratio as a function of position
        """
        #get variables, based on daterange
        if daterange:
            peak_ratio = ( self.matches_date_range['peak_flux_ap'] /
                           self.matches_date_range['int_flux_nvss'] )
            delta_ra = self.matches_date_range['delta_ra']
            delta_dec = self.matches_date_range['delta_dec']
            figname = ("{0}_B{1}_{2}_{3}_"
                       "peak_ratio_position.png").format(self.pbname,self.beam,
                                                         self.startdate,self.enddate)
        else:
            peak_ratio = self.matches['peak_flux_ap']/self.matches['int_flux_nvss']
            delta_ra = self.matches['delta_ra']
            delta_dec = self.matches['delta_dec']
            figname = "{0}_B{1}_peak_ratio_position.png".format(self.pbname,self.beam)
        #setup figure
        fig, ax = plt.subplots(figsize=(8,8))
        #scatter plot w/ color mapped to peak_ratio
        #later add size proportional to source size?
        sc= ax.scatter(delta_ra,delta_dec,
                       c=peak_ratio,s=20)
        #set color bar
        fig.colorbar(sc,ax=ax)
        #set titles and labels
        ax.set_title('Apertif peak flux / NVSS integrated flux')
        ax.set_xlabel('Delta RA [arcmin]')
        ax.set_ylable('Delta Dec [arcmin]')
        #save plot to pre-defined output, based on PB/CB
        figpath = os.path.join(figdir,figname)
        plt.savefig(figpath)

    def int_ratio_position(self,daterange=False):
        """
        Plot ratio of integrated flux as a function of position
        """
        #get variables, based on date range
        if daterange:
            int_ratio = ( self.matches_date_range['int_flux_ap'] /
                           self.matches_date_range['int_flux_nvss'] )
            delta_ra = self.matches_date_range['delta_ra']
            delta_dec = self.matches_date_range['delta_dec']
            figname = ("{0}_B{1}_{2}_{3}_"
                       "int_ratio_position.png").format(self.pbname,self.beam,
                                                         self.startdate,self.enddate)
        else:
            int_ratio = self.matches['int_flux_ap']/self.matches['int_flux_nvss']
            delta_ra = self.matches['delta_ra']
            delta_dec = self.matches['delta_dec']
            figname = "{0}_B{1}_int_ratio_position.png".format(self.pbname,self.beam)
        #setup figure
        fig, ax = plt.subplots(figsize=(8,8))
        #scatter plot w/ color mapped to peak_ratio
        #later add size proportional to source size?
        sc = ax.scatter(delta_ra,delta_dec,
                        c=int_ratio,s=20)
        #set color bar
        fig.colorbar(sc,ax=ax)
        #set titles and labels
        ax.set_title('Apertif integrated flux / NVSS integrated flux')
        ax.set_xlabel('Delta RA [arcmin]')
        ax.set_ylable('Delta Dec [arcmin]')
        #save plot to pre-defined output, based on PB/CB
        figpath = os.path.join(figdir,figname)
        plt.savefig(figpath)
