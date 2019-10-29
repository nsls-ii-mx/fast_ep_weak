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
import argparse

from xml.dom import minidom

parser = argparse.ArgumentParser()
parser.add_argument("a")
args = parser.parse_args()

print ("Your data is in: "), args.a

os.system("fast_dp %s -j 20 -1 10 -N 1100" % args.a)

# Parsing fast_dp.xml
if(os.path.isfile('fast_dp.xml')):
  xmldoc_fast_dp = ET.parse("fast_dp.xml")
  root_fast_dp = xmldoc_fast_dp.getroot()
  root0c = root_fast_dp[0].getchildren()
  filePath   = root0c[1].find("filePath").text
  #filePath   = root_fast_dp[0][1].find("filePath").text
  print ("Your files are in: "), filePath
  spacegroup = root_fast_dp[1].find("spaceGroup").text
  print ("Space group is: "), spacegroup
  rMerge = root_fast_dp[2][1].find("rMerge").text
  print ("R merge is: "), rMerge
  ResLo = root_fast_dp[2][1].find("resolutionLimitLow").text
  ResHi = root_fast_dp[2][1].find("resolutionLimitHigh").text
  print ("Resolution limits:\n   Low:  "), ResLo 
  print ("   High: "), ResHi 

percent_solvent = 0.50
asu_mol = 1
solvent = 0.50

# Parsing MATTHEWS_COEF.xml
if(os.path.isfile('MATTHEWS_COEF.xml')):
  xmldoc_matthews = ET.parse("MATTHEWS_COEF.xml")
  root_matthews = xmldoc_matthews.getroot()
  for result_all in root_matthews.findall("result"):
     asu_mol = float(result_all.get("nmol_in_asu"))
     print ("ASU molecules: "), asu_mol
     solvent = float(result_all.get("percent_solvent"))/100
     print ("Percent solvent: "), solvent

#Writing xml file

#import xml.etree.cElementTree as ET

root = ET.Element("root")
CompInData = ET.SubElement(root, "CompInData")
OutputData = ET.SubElement(root, "OutputData")

ET.SubElement(CompInData, "SolventContent").text = str(solvent)

ET.SubElement(OutputData, "SpaceGroup").text = str(spacegroup)
ET.SubElement(OutputData, "Rmerge").text = str(rMerge)
ET.SubElement(OutputData, "ResolutionLimitLo").text = str(ResLo)
ET.SubElement(OutputData, "ResolutionLimitHi").text = str(ResHi)
ET.SubElement(OutputData, "MoleculesPerASU").text = str(asu_mol)

tree = ET.ElementTree(root)
tree.write("pipelineoutput.xml")
