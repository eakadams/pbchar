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
                 matchfile):
        """
        Initialize PB object
        Inputs:
        - matchfile (str): Path of matchfile with all cross-match info
        """
        #read in match file
        self.matches = ascii.read(matchfile,format='csv')
        basename = os.path.basename(matchfile)
        names = basename.split("_")
        pbname = names[0]
        names2 = names[1].split(".")
        beam = names2[0]
        print(pbname,beam)

    def go(self):
        """
        Run all available analysis plots:
        - peak_ratio_position
        - int_ratio_position

        Optionally, same plots also date-range-limited:
        - peak_ratio_position_date_range
        - int_ratio_position
        """
        #have to write functions for plots first
        pass


    def peak_ratio_position(self):
        """
        Plot peak ratio as a function of position
        """
        #setup figure
        fig, ax = plt.subplots(figsize=(8,8))
        #get peak ratio
        peak_ratio = self.matches['peak_flux_ap']/self.matches['int_flux_nvss']
        #scatter plot w/ color mapped to peak_ratio
        #later add size proportional to source size?
        ax.scatter(self.matches['delta_ra'],self.matches['delta_dec'],
                   c=peak_ratio,s=20)
        #save plot to pre-defined output, based on PB/CB
        figname = "{0}_B{1}_peak_ratio_position.png".format(self.pbname,self.beam)
        figpath = os.path.join(figdir,figname)
        plt.savefig(figpath)
