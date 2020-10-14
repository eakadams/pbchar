#PB char for single beam

from __future__ import print_function

__author__="E.A.K. Adams"

"""
Code for doing PB correction of a single beam
and comparing to NVSS, writing results to a 
text file. Create a Beam object that handles
everything. 

Improvement to-do list:
- Make NVSS filter size a param to init Beam
- helper function for pbfits path (handle different types)
"""

import os
import subprocess
import shutil
import bdsf
from astroquery.vizier import Vizier
from apercal.libs import lib
from astropy.coordinates import SkyCoord
from astropy.io import ascii
from astropy.io import fits
from astropy.table import Table
from astropy import units as u


#import glob
#import sys
#import numpy as np





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
                 skipcheck = False,
                 masklevel = None):
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
        - masklevel: Level at which to mask image after PB correction; 
                     e.g. 0.2 for 20% of PB response
                     default is None
        """
        #add initial values
        self.taskid = str(taskid)
        self.beam = '{:02d}'.format(beam)
        self.pbdir = pbdir
        self.pbname = pbname
        self.skipcheck = skipcheck
        #check that masklevel is a float
        if type(masklevel) is float:
            self.masklevel = masklevel
        else:
            self.masklevel = None
            print("Mask level not a float; no masking will be done")

        #setup output locations
        #should include info about PBs being used and beam
        #taskids go in file names
        self.outputdir = os.path.join(outputdir,self.pbname,self.beam)
        #make outputdir if it doesn't exist
        if not os.path.isdir(self.outputdir):
            os.makedirs(self.outputdir)

        #working dir
        #needs to have taskid/beam in case I run things in parallel
        #want to be able to clean up working dir at end and not
        #worry about writing conflicts
        self.workingdir = os.path.join(workingdir,self.taskid,self.beam)
        #make working dir if it doesn't exist
        if not os.path.isdir(self.workingdir):
            os.makedirs(self.workingdir)

        #define paths I will want later
        self.fitspath = None #can only define this when know image name
        #m for "map"
        self.impath = os.path.join(self.workingdir,"m")
        #sm for "SMoothed Map"
        self.smimpath = os.path.join(self.workingdir,"sm")
        #pbr for "Primary Beam Regridded
        self.pbpath = os.path.join(self.workingdir,'pbr')
        #sm_pb for "SMoothed and PB corrected"
        self.pbsmimpath = os.path.join(self.workingdir,"sm_pb")
        
        # pbsmfits fits; sm_pb_<taskid>_mfs.fits in output dir; may want to keep
        #but not for now while testing
        self.pbsmfits = os.path.join(self.outputdir,"sm_pb_{}_mfs.fits".format(self.taskid))
        self.pbsmfits = os.path.join(self.workingdir,"sm_pb.fits")
        #soruce finding results; also to output dir, may want to keep
        #but not for now while testing
        self.bdsf_file = os.path.join(self.workingdir,"bdsf.sav")
        self.bdsf_output = os.path.join(self.workingdir,"bdsf_output.srl")

        #final output; csv file with matches
        self.match_output = os.path.join(self.outputdir,"{0}_matches.csv".format(self.taskid))

        #variables for RA, Dec of center of image
        self.ra = None
        self.dec = None

        #nvss table
        self.nvss_table = None
        #cross-match table
        self.match_table = None
        
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
        else:
            print('Output already exists; not running code')
            
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
        elif os.path.exists(self.match_output):
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
                print(self.fitspath)
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
        #change to check fitspath exist
        #setting status but not really using now
        if os.path.exists(self.fitspath):
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
                    print(e)

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
        #check that image exists
        if os.path.isdir(self.impath):    
            #check that output doesn't exist
            if not os.path.isdir(self.smimpath):
                #do convolution
                convol = lib.miriad('convol')
                convol.map = self.impath
                convol.out = self.smimpath
                convol.fwhm = '45,45' #NVSS resolution
                convol.pa = '0' #should be default but had case where not in fits header
                convol.options = 'final'
                try:
                    convol.go()
                except Exception as e:
                    #print(e)
                    self.status = False
                    print(("Convoling beam {0} of taskid {1} failed").
                          format(self.beam,self.taskid))
                    print(e)
        else:
            #image to convolve doesn't exist
            self.status = False
            self.smimpath = None
        

    def regrid_pb(self):
        """
        Convert PB image to miriad in working directory
        Regrid image
        """
        #first get PB image and write to miriad in workingdir
        #naming convention might depend on type of PB image
        #drift, Alexander's method, etc
        #check if drift
        if len(self.pbname) == 6:
            #YYMMDD name used for driftscans
            pbfits = os.path.join(self.pbdir,"{0}_{1}_I_model_reg.fits".
                                  format(self.pbname,self.beam))
            #may have to add a check about string length, but wait until it's a problem

        #having defined pbfits above,
        #check it exists, retrieve to miriad in workingdir
        #also check that output doesn't exist
        pbim = os.path.join(self.workingdir,"pb")
        if (os.path.exists(pbfits) and not os.path.isdir(pbim)):
            #get pbimage to miriad           
            fits = lib.miriad('fits')
            fits.op = 'xyin'
            fits.in_ = pbfits
            fits.out = pbim
            try:
                fits.go()
            except Exception as e:
                self.status = False
                print(("Conversion to miriad image failed "
                       "for primary beam {0}").
                          format(self.beam))
                print(e)

        #this seems to then have regrid do nothing
        #but it lets me run code and can update/fix this later
        #get projection center of continuum image for regridding
        #and update center of PB image manually
        #is this the right thing to do?
        #check w/ Thijs & Helga.....
        if (os.path.isdir(self.smimpath)
            and
            os.path.isdir(pbim)):
            #get center values
            gethd = lib.miriad('gethd')
            gethd.in_ = os.path.join(self.smimpath,'crval1')
            #gethd.format = 'hms'
            ra_ref = gethd.go()
            gethd.in_ = os.path.join(self.smimpath,'crval2')
            #gethd.format = 'dms'
            dec_ref = gethd.go()
            #update center
            puthd = lib.miriad('puthd')
            puthd.in_ = os.path.join(pbim,'crval1')
            puthd.value = ra_ref[0]
            puthd.type = 'double'
            puthd.go()
            puthd.in_ = os.path.join(pbim,'crval2')
            puthd.value = dec_ref[0]
            puthd.type = 'double'
            puthd.go()
            #keep ra & dec in object; will need later
            self.ra = ra_ref[0] * u.rad
            self.dec = dec_ref[0] * u.rad
        else:
            #getting projection center of (smoothed) image failed
            self.status = False
            print(("Updating projection center from beam {0}, "
                   "based on taskid {1}, failed").
                  format(self.beam,self.taskid))

        #regrid image
        #define out put, pbr - PB R egridded
        
        #check pbim exists for regridding
        #plus smoothed path for template
        #also that output doesn't exist
        if (os.path.isdir(pbim) and
            os.path.isdir(self.smimpath) and
            not os.path.isdir(self.pbpath)):
            regrid = lib.miriad('regrid')
            regrid.in_ = pbim
            regrid.out = self.pbpath
            regrid.tin = self.smimpath
            regrid.axes = '1,2'
            #what if I try options = offset
            #then maybe keep image values but on right grid
            regrid.options = 'offset'
            #then my image is resize so things don't work later
            #what if I add desc value?
            #use helper function to get values
            ###ra_vals = self.get_hdr('ra',self.smimpath)
            ###dec_vals = self.get_hdr('dec',self.smimpath)
            ###desc = "{0},{1},{2},{3},{4},{5},{6},{7}".format(
            ###    ra_vals[0],ra_vals[1],ra_vals[2],ra_vals[3],
            ###    dec_vals[0],dec_vals[1],dec_vals[2],dec_vals[3])
            try:
                regrid.go()
            except Exception as e:
                self.status = False
                print(("Regridding failed for primary beam {0},"
                       "reference taskid {1}").format(self.beam,self.taskid))
                print(e)
        else:
            self.status = False
            print(("Either primary beam or smoothed continuum "
                   "image missing for beam {0}, taskid{1}").
                  format(self.beam,self.taskid))
            

    def get_hdr(self,axis,image):
        """
        Get four sets of header values for given axis from image
        CRVAl, CRPIX, CDELT, NAXIS
        """
        if axis is 'ra':
            num = '1'
        elif axis is 'dec':
            num = '2'
        else:
            print("axis not recognized; default first axis")
            num = '1'

        gethd = lib.miriad('gethd')

        gethd.in_ = os.path.join(image,'crval{}'.format(num))
        crval = gethd.go()

        gethd.in_ = os.path.join(image,'crpix{}'.format(num))
        crpix = gethd.go()

        gethd.in_ = os.path.join(image,'cdelt{}'.format(num))
        cdelt = gethd.go()

        gethd.in_ = os.path.join(image,'naxis{}'.format(num))
        naxis = gethd.go()
            
        return crval,crpix,cdelt,naxis

    def do_pb(self):
        """
        Apply primary beam correction using immath
        """
        #do PB correction
        #check that regridded PB image and smoothed cont image exist
        #also that output doesn't
        if (os.path.isdir(self.smimpath) and
            os.path.isdir(self.pbpath) and
            not os.path.isdir(self.pbsmimpath)):
            maths = lib.miriad('maths')
            maths.exp = "'<{0}>/<{1}>'".format(self.smimpath, self.pbpath)
            maths.out = self.pbsmimpath
            print(maths.exp)
            print(maths.out)
            #add a mask at set level of PB response
            if self.masklevel is not None:
                #add a mask
                maths.mask = "'<{0}>.gt.{1}'".format(self.pbpath,self.masklevel)
                print(maths.mask)
            try:
                maths.go()
            except Exception as e:
                self.status = False
                print(("Primary beam correction failed for "
                       "beam {0}, taskid {1}").format(self.beam,self.taskid))
                print(e)

        if (os.path.isdir(self.pbsmimpath) and
            not os.path.exists(self.pbsmfits)):
            #also write to fits file; needed for pyBDSF
            fits = lib.miriad('fits')
            fits.op = 'xyout'
            fits.in_ = self.pbsmimpath
            fits.out = self.pbsmfits
            try:
                fits.go()
            except Exception as e:
                self.status = False
                print(("Writing cont image to fits for source finding failed "
                       "for beam {0} of taskid {1}").format(self.beam,self.taskid))
                print(e)
                
            

    def find_sources(self):
        """
        Find sources in smoothed & pb-corrected image
        Write them to a temp file in working directory
        """
        #dont' do source finding if output already exists
        if not os.path.exists(self.bdsf_output) and not self.skipcheck:
            bdsf.process_image(self.pbsmfits,adaptive_rms_box=True,
                               thresh_isl=3.0, thresh_pix=5.0).write_catalog(outfile=self.bdsf_output,
                                                                             format='ascii',
                                                                             clobber=True)
            
        
        #add something like self.sourcepath
        #if I can keep sources in a Table add add that to ibject
        #that would also be great

    def get_nvss(self):
        """
        Get NVSS sources
        Using astroquery Vizier functionality
        """
        #NVSS catalog identifier
        nvsscat = 'VIII/65/nvss'
        #get rid of row limit
        Vizier.ROW_LIMIT = -1
        #get sky coord of cont image
        v = Vizier(columns = ['NVSS','RAJ2000','DEJ2000','S1.4','e_S1.4',
                              'MajAxis','MinAxis'],
                   row_limit = -1,
                   column_filters = {"MajAxis":"<50", "MinAxis":"<50"})
        print(self.ra.to(u.deg),self.dec.to(u.deg))
        result = v.query_region(SkyCoord(self.ra.to(u.deg),
                                         self.dec.to(u.deg),
                                         frame='icrs',
                                         unit = (u.deg,u.deg)),
                                radius=60*u.arcmin,
                                catalog=nvsscat)
        #print(result)
        #print(result[nvsscat])

        self.nvss_table = result[nvsscat]
        #print(self.nvss_table.colnames)

        #print(self.nvss_table['MajAxis','MinAxis'])
        #add something like self.nvss_sources which is table

    def cross_match(self):
        """
        Do the cross-matching
        Filtered on size in retrieval above
        Make sure coordiante are close
        """
        #read in source list from pybdsf output
        bdsf_sources = ascii.read(self.bdsf_output,guess=True,
                                  data_start=5,header_start=4)
        print(bdsf_sources.colnames)

        #initialize empty lists to hold output from successful cross-matches
        peak_flux_ap = [] #apertif
        int_flux_ap = []
        int_flux_nvss = []
        deltara = [] #delta RA (Apertif) from image center
        deltadec = [] #delta dec (Apertif) from image center
        radius = [] #distance from field center
        pb_level = [] #primary beam response at same position

        #get NVSS skycoord object
        #print(self.nvss_table['RAJ2000','DEJ2000'][0:10])
        nvss_coords = SkyCoord(ra=self.nvss_table['RAJ2000'],
                               dec=self.nvss_table['DEJ2000'],
                               frame='icrs',
                               unit=(u.hourangle,u.deg))

        #print(nvss_coords[0:10])
        #get center  iamge coord also
        center_coord  = SkyCoord(ra=self.ra.to(u.deg),
                                 dec=self.dec.to(u.deg),
                                 unit=(u.deg,u.deg),
                                 frame='icrs')

        #open primary beam image for getting pb level
        #might want a helper function for getting pbfits in future....
        if len(self.pbname) == 6:
            #YYMMDD name used for driftscans
            pbfits = os.path.join(self.pbdir,"{0}_{1}_I_model_reg.fits".
                                  format(self.pbname,self.beam))

        with fits.open(pbfits) as hdul:
            pbdata = hdul[0].data

            
        #iterate through every Apertif sources
        for i in range(len(bdsf_sources)):
            #get skycoord of source
            source_coord = SkyCoord(ra=bdsf_sources['RA'][i],dec=bdsf_sources['DEC'][i],
                                    unit=(u.deg,u.deg))
            #get closest match in NVSS catalog
            idx, sep2d, dist3d = source_coord.match_to_catalog_sky(nvss_coords)
            #print(sep2d.to(u.arcsec).value)
            #check if separation is w/in 5"
            if sep2d  < 5*u.arcsec:
                #append values to list
                peak_flux_ap.append(bdsf_sources['Peak_flux'][i])
                int_flux_ap.append(bdsf_sources['Total_flux'][i])
                int_flux_nvss.append((self.nvss_table['S1.4'][idx])/1000.) #record in  Jy, match bdsf
                d_ra, d_dec = center_coord.spherical_offsets_to(source_coord)
                r = center_coord.separation(source_coord)
                deltara.append(d_ra.to(u.arcmin).value)
                deltadec.append(d_dec.to(u.arcmin).value)
                radius.append(r.to(u.arcmin).value)
                xpix = int(bdsf_sources['Xposn'][i]) - 1 #0-index; force integer pixel
                ypix = int(bdsf_sources['Yposn'][i]) - 1 #0-index; force interger pixel
                pbval = pbdata[ypix,xpix] #axes reversed
                pb_level.append(pbval)

        #print(flux_ap,flux_nvss,radius)

        #next create a table that is cross match
        print(peak_flux_ap)
        self.match_table = Table([peak_flux_ap, int_flux_ap,
                                  int_flux_nvss,
                                  deltara,deltadec,radius,pb_level],
                                 names=('peak_flux_ap','int_flux_ap',
                                        'int_flux_nvss','delta_ra',
                                        'delta_dec','radius','pb_level'))

            

    def record(self):
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
        #write out file
        ascii.write(self.match_table,self.match_output,format='csv',overwrite=True)
        #self.match_table.write(self.match_output,overwrite=True,format='csv')

        #write out files

    def clean_up(self):
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
        #remove working directory to remove all intermediate files
        shutil.rmtree(self.workingdir,ignore_errors=True)
