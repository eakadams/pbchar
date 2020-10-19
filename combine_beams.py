#python code to combine beam information in one master file

from __future__ import print_function

"""
Combine information from all taskids
in one master file for each  primary beam
Write these to pbchar directory so that they
can be part of git repository and accessed 
in other locations for analysis that comes later
"""

from astropy.io import ascii
from multiprocessing import Pool
import argparse
import os
import glob
from astropy.table import Table

#get global level paths
this_dir,this_filename = os.path.split(__file__)
pbchardir = this_dir
filedir = os.path.join(pbchardir,"files")

parser  = argparse.ArgumentParser(description='Collect information per beam for given PB set')
parser.add_argumentt("PBdir", help='PB directory w/ cross-match info')

def collect_info(beam):
    """
    Function to collect information for a single beam
    """
    beamdir = os.path.join(PBdir,"{0:02d}".format(beam))
    #get list of csv files for each tid
    tid_csv_list = glob.glob(os.path.joing(beamdir,"*matches.csv"))
    #set up empty arrays to hold things
    tid_array = np.array([])
    peak_flux_array = np.array([])
    int_flux_array = np.array([])
    int_flux_nvss_array = np.array([])
    delta_ra_array = np.array([])
    delta_dec_array = np.array([])
    radius_array = np.array([])
    pb_level_array = np.array([])
    #iterate through matches and append to arrays
    for tid_csv in tid_csv_list:
        tid_csv_name = os.path.basename(tid_csv)
        tid = tid_csv[0:9]
        t = ascii.read(tid_csv,format='csv')
        tmp_tid_array = np.full(len(t),tid)
        tid_array = np.append(tid_array,tmp_tid_array)
        peak_flux_array = np.append(peak_flux_array,t['peak_flux_ap'])
        int_flux_array = np.append(int_flux_array,t['int_flux_ap'])
        int_flux_nvss_array = np.append(int_flux_nvss_array,t['int_flux_nvss'])
        delta_ra_array = np.append(delta_ra_array,t['delta_ra'])
        delta_dec_array = np.append(delta_dec_array,t['delta_dec'])
        radius_array = np.append(radius_array,t['radius'])
        pb_level_array = np.append(pb_level_array,t['pb_level'])
    #put info into a table
    combined_beam_table = Table([tid_array,peak_flux_array,int_flux_array,
                                int_flux_nvss_array,delta_ra_array,
                                delta_dec_array,radius_array,pb_level_array],
                                names = ('ObsID','peak_flux_ap','int_flux_ap',
                                            'int_flux_nvss','delta_ra',
                                            'delta_dec','radius','pb_level'))
    #write table to pbchar, filedir
    pbname = os.path.basename(PBdir)
    output = os.path.join(filedir,"{0}_{1:02d}.csv".format(pbname,beam))
    ascii.write(combined_beam_table,output,format='csv',overwrite=True)
    
    
if __name__ == '__main__':
    #do parallelization, per beam
    beams = np.arange(40)
    #set up pool; 10 cores
    pool = Pool(10)
    #run by mapping onto pool
    pool.map(collect_info,beams)
    #close and join
    pool.close()
    pool.join()
