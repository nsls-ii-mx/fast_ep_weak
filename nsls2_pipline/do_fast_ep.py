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

print ("Your mtz file is: "), args.a
print ("fast_ep sad=%s" % args.a)

#if(os.path.isfile("fast_dp.mtz")):

os.system("fast_ep sad=%s" % args.a)
os.system("./do_fast_ep_2_xml.com")

eptree = ET.parse('fast_ep.xml')
best_hand = eptree.findall('./OutputData/BestHand')[0].text
n_sites = eptree.findall('./OutputData/NSites')[0].text

print ("Best hand is:       "), best_hand
print ("Number of sites is: "), n_sites

#Writing xml file

tree = ET.parse('pipelineoutput.xml')
root = tree.getroot()

OutputData = root.find('OutputData')
CompInData = root.find('CompInData')
ET.SubElement(OutputData, "BestHand").text = str(best_hand)
ET.SubElement(CompInData, "NSites").text = str(n_sites)

print ET.tostring(root)

tree.write("pipelineoutput.xml")
