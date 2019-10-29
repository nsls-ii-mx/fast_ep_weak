import sys
import os
import time
import copy
import exceptions
import traceback
import subprocess
import os.path
import xml.etree.ElementTree as ET
import xmltodict
import math
import logging
from xml.dom import minidom

data_input = raw_input("Enter fast_dp data location: ")

#Run fast_dp
if(os.path.isfile(data_input)):
  print ("\nRunning fast_dp on data: %s" % data_input)
  os.system("python do_fast_dp.py %s" % data_input)
else:
  print "FAILED to find raw input data"

#Run fast_ep
#Manual:  python do_fast_ep.py fast_dp.mtz
if(os.path.isfile("fast_dp.mtz")):
  print ("Running fast_ep on data: fast_dp.mtz")
  os.system("python do_fast_ep.py fast_dp.mtz")
else:
  print "FAILED to find fast_dp output"

#Read output files
if(os.path.isfile('sad.mtz')):
  print ("Parsing data using: sad.mtz")
  os.system("python do_sad_parse.py sad.mtz")
else:
  print "FAILED to find fast_ep output"

#Read solve structure
if(os.path.isfile('phs.mtz')):
  print ("Parsing data using: sad.mtz")
  os.system("python final.py")
else:
  print "FAILED to find SHELX output"
