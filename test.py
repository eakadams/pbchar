#command line code for testing

from __future__ import print_function

import beam as beam

print("Load Beam object; 190807041, 1; mask 10%")
B = beam.Beam(190807041,1,masklevel=0.1,workingdir='tmp')

print("Load Beam object; 190823041, 9; mask 10%")
B = beam.Beam(190823041,9,masklevel=0.1,workingdir='tmp')

print("Get cont image")
B.get_cont_image()

print("Smooth cont image")
B.convolve()

print("Get and regrid PB")
B.regrid_pb()

print("Do PB correction")
B.do_pb()

print("Test done")
