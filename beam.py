#PB char for single beam

from __future__ import print_function

__author__="E.A.K. Adams"

"""
Code for doing PB correction of a single beam
and comparing to NVSS, writing results to a 
text file. Create a Beam object that handles
everything. 
"""

import os
import subprocess

from apercal.libs import lib

import glob
import sys
from astropy.table import Table
import numpy as np
from astropy.io import ascii
from astropy.coordinates import SkyCoord
from astropy import units as u



class Beam(object):

    #stole from mosaic code, used for subprocess call checking
    FNULL = open(os.devnull, 'w')
    
    def __init__(self,
                 taskid,
                 beam,
                 pbdir = '/tank/apertif/driftscans/DR1/190912/continuum',
                 pbname = '190912',
                 outputdir = '/tank/adams/pbchar/',
                 workingdir = '/tank/adams/pbchar/tmp/',
                 skipcheck = False):
        """
        Initialize Beam object
        Inputs:
        - taskid : Taskid of observation
        - beam: Beam number of observation
        - pbdir: Path to primary beam images; default is DR1 on happili-05
        - pbname (string) : Name describing PB images used
        - outputdir: Path for location to write output; 
                     default is personal space on happili-05
        - workingdir: Path for writing files as part of process
        - skipcheck: Can skip checking for output to force a rerun
        """
        #add initial values
        self.taskid = str(taskid)
        self.beam = '{:02d}'.format(beam)
        self.pbdir = pbdir
        self.pbname = pbname

        #setup output locations
        #should include info about PBs being used and beam
        #taskids go in file names
        self.outputdir = os.path.join(outputdir,self.pbname,self.beam)

        #working dir
        self.workingdir = workingdir

        #setup a status
        #can query to see if steps should be run
        #also useful output if things were successful overall
        self.status = True
        
    def go(self):
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
            print('Doing PB correction and source finding')
            self.get_cont_image()
            self.convolve()
            self.regrid_pb()
            self.do_pb()
            self.find_sources()
            self.get_nvss()
            self.cross_match()
            self.record()
            self.clean_up()
            

    def check(self):
        """
        Do a check to see if beam has already been run
        Then don't repeat
        Can also force the check to be skipped if want to
        redo things on purpose
        """
        if self.skipcheck:
            check = False
        #do the check
        #simple - does output directory exist?
        elif os.path.exists(self.outputdir):
            check = True
        else:
            check = False

        return check


    def get_cont_image(self):
        """
        Get continuum image specified by taskid, beam
        Want to retrieve from ALTA
        Can use mosaic module in Apercal as a guide to how to do this
        Convert to miriad format
        """
        #directory for beam/taskid
        altadir = ("/altaZone/archive/apertif_main/"
                   "visibilities_default/{0}_AP_B{1}").format(
                        self.taskid, self.beam.zfill(3))

        #check if beam/taskid available in ALTA
        if self.check_alta_path(altadir) == 0:
            #found beam so can continue
            #have to find continuum image name since it can be different
            fits_name = ''
            alta_fits_path = ''

            for k in range(10):
                fits_name = "image_mf_{0:02d}.fits".format(k)
                alta_fits_path = os.path.join(altadir, fits_name)
                if self.check_alta_path(alta_fits_path) == 0:
                    break
                else:
                    # make empty again when no image was found
                    fits_name = ''
                    continue
            if fits_name == '':
                #failed to find image. Note this and print warning
                self.status = False
                ("No continuum image found on ALTA for "
                 "beam {0} of taskid {1}").format(self.beam, self.taskid)
            else:
                #create working directory
                if not os.path.exists(self.workingdir):
                    os.makedirs(self.workingdir)
                #check if file is already there
                self.fitspath = os.path.join(self.workingdir,fits_name)
                if not os.path.exists(self.fitspath):
                    #copy continuum image to this directory
                    return_msg = self.get_data_from_alta(alta_fits_path,
                                                        self.fitspath)
                    if return_msg != 0.0:
                        print(("Failed to retrieve continuum "
                               "image for beam {0} of taskid {1}").
                              format(self.beam,self.taskid))
                        self.status = False

        #now, as long as status is good, convert to miriad
        if self.status is True:
            #set miriad image path
            #keep it short because of string length issues
            #means will have to be good about doing cleanup!
            self.impath = os.path.join(self.workingdir,
                                       ("image_{0}.map").
                                       format(self.beam))
            #check that miriad image doesn't already exist and do conversion
            if not os.path.isdir(self.impath):
                fits = lib.miriad('fits')
                fits.op = 'xyin'
                fits.in_ = self.fitspath
                fits.out = self.impath
                try:
                    fits.go()
                except Exception as e:
                    self.status = False
                    print(("Conversion to miriad image failed "
                           "for beam {0} of taskid {1}").
                          format(self.beam,self.taskid))

    def check_alta_path(self,alta_path):
        """
        Need to test if ALTA path exists
        Code from mosaic_v2 module in Apercal (by R. Schulz, prob)
        """
        alta_cmd = "ils {}".format(alta_path)
        return_msg = subprocess.call(alta_cmd, shell=True,
                                     stdout=self.FNULL, stderr=self.FNULL)
        return return_msg

    def get_data_from_alta(self, alta_file_name, output_path):
        """
        Function to get files from ALTA
        Code from mosaic_v2 module in Apercal
        """
        # set the irod files location
        irods_status_file = os.path.join(
            os.getcwd(),
            "transfer_{}_img-icat.irods-status".format(
                os.path.basename(alta_file_name).split(".")[0]))
        irods_status_lf_file = os.path.join(
            os.getcwd(),
            "transfer_{}_img-icat.lf-irods-status".format(
                os.path.basename(alta_file_name).split(".")[0]))

        # get the file from alta
        alta_cmd = ("iget -rfPIT -X {0} --lfrestart {1} "
                    "--retries 5 {2} {3}").format(
                        irods_status_file, irods_status_lf_file,
                        alta_file_name, output_path)
        return_msg = subprocess.check_call(
            alta_cmd, shell=True, stdout=self.FNULL, stderr=self.FNULL)

        return return_msg


    def convolve(self):
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
