#PB char for single beam

from __future__ import print_function

__author__="E.A.K. Adams"

"""
Code for doing PB correction of a single beam
and comparing to NVSS, writing results to a 
text file. Create a Beam object that handles
everything. 
"""

import glob
import sys
from astropy.table import Table
import numpy as np
from astropy.io import ascii
import os
import shutil
import matplotlib.pyplot as plt
from astropy.coordinates import SkyCoord
from astropy import units as u
from apercal.libs import lib


class Beam(object):
    def __init__(self,
                 taskid,
                 beam,
                 pbdir = '/tank/apertif/driftscans/DR1/190912/continuum'
                 outputdir = '/tank/adams/pbchar/',
                 workingdir = '/tank/adams/pbchar/tmp/',
                 skipcheck = False):
        """
        Initialize Beam object
        Inputs:
        - taskid : Taskid of observation
        - beam: Beam number of observation
        - pbdir: Path to primary beam images; default is DR1 on happili-05
        - outputdir: Path for location to write output; 
                     default is personal space on happili-05
        - workingdir: Path for writing files as part of process
        - skipcheck: Can skip checking for output to force a rerun
        """
        #add initial values
        self.taskid = taskid
        self.beam = beam
        self.pbdir = pbdir

        #setup output locations
        #should include info about PBs being used....
        self.outputdir = outputdir

    def go():
        """
        Run the body of code:
        (0. Check if information already exists)
        1. Retrieve continuum image
        2. Convolve to NVSS resoltuion
        4. Regrid PB image to cont image
        4. Apply PB correction
        5. Run source finder
        6. Get NVSS sources over same region
        7. Crossmatch sources
        8. Record crossmatches and fluxes
        9. Clean up
        """

        #check if files already exist
        #if so, don't need to redo things
        checkresult = self.check()

        if not checkresult:
            self.get_cont_image()
            self.convolve()
            self.regrid_pb()
            self.do_pb()
            self.find_sources()
            self.get_nvss()
            self.cross_match()
            self.record()
            self.clean_up()
            

    def check():
        """
        Do a check to see if beam has already been run
        Then don't repeat
        Can also force the check to be skipped if want to
        redo things on purpose
        """
        if self.skipcheck:
            check = False
        else:
            #do the check!
            check = False

        return check


    def get_cont_image():
        """
        Get continuum image specified by taskid, beam
        Want to retrieve from ALTA
        Can use mosaic module in Apercal as a guide to how to do this
        Convert to miriad format
        """
        altadir = ("/altaZone/archive/apertif_main/"
                   "visibilities_default/{0}_AP_B{1}").format(
                        self.taskid, self.beam.zfill(3))
        
        
        #add something like self.impath

    def convolve():
        """
        Convolve to NVSS resolution
        """
        #do convolution

        #add something like self.smimpath
        

    def regrid_pb():
        """
        Convert PB image to miriad in working directory
        Regrid image
        """

        #check if working directory exists, make if needed

        #regrid image

        #add something like self.pbbath


    def do_pb():
        """
        Apply primary beam correction using immath
        """
        #do PB correction

        #add something like self.pbsmimpath


    def find_sources():
        """
        Find sources in smoothed & pb-corrected image
        Write them to a temp file in working directory
        """

        #add something like self.sourcepath
        #if I can keep sources in a Table add add that to ibject
        #that would also be great

    def get_nvss():
        """
        Get NVSS sources
        """

        #have to decide whether it's better to query
        #or have catalog locally
        #probablyb etter to  query because then it's smaller
        #for source matching
        #and I like idea of using tools for querying

        #add something like self.nvss_sources which is table

    def cross_match():
        """
        Do the cross-matching
        Apply filters based on:
        - flux?
        - size?
        - flux + PB location?
        """

        #add something like self.matches which is a table of matches

    def record():
        """
        Write out matches to text files for later access
        Format is something like:
        outputdir / beam / PB+taskid
        Want to include PB date info to allow comparisons between 
        different  PB images
        Could alternatively think  of that as a different directory:
        outputdir / PB / beam / taskid.csv
        That might be a better, cleaner organization
        """

        #write out files

    def clean_up():
        """
        Clean up any temporary files
        This will be:
        - miriad PB
        - regridded PB
        - retrieved cont image (fits & miriad)
        - convolved cont image
        - PB-corrected cont image
        - pybdsf output (if text files)
        - NVSS catalog (if text files)
        Do I want to optionally keep smoothed & pb-corrected continuum image?
        Could allow starting from that point and just doing later steps
        Would have to modify starting check
        But could allow updates if I improve later steps
        Not sure about total size, but probably small
        I would want to write to fits to keep though, so seems complicted
        Maybe future addition?
        """
