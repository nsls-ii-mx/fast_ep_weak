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

print ("Parsing data set: "), args.a

# Parsing fast_dp.xml
xmldoc_fast_dp = ET.parse("fast_dp.xml")
root = xmldoc_fast_dp.getroot()
spacegroup = root[1].find("spaceGroup").text
SG = spacegroup.replace(" ", "") # removes the white spaces
a_cell = root[1].find("refinedCell_a").text
b_cell = root[1].find("refinedCell_b").text
c_cell = root[1].find("refinedCell_c").text
alpha = root[1].find("refinedCell_alpha").text
beta = root[1].find("refinedCell_beta").text
gamma = root[1].find("refinedCell_gamma").text
print ("Space group: "), SG
print ("Cell constants:\n   a      "), a_cell
print ("   b      "), b_cell
print ("   c      "), c_cell
print ("   alpha  "), alpha
print ("   beta   "), beta
print ("   gamma  "), gamma

# Defining fast_ep.log and generating sh input files for shelxe
# and dm.
fast_ep_log = open("fast_ep.log", "r")
fast_ep_lines = fast_ep_log.readlines()
fast_ep_string = str(fast_ep_lines)

nsitesFastep = 10
# Parse nsites from fast_ep.log
for line in open("fast_ep.log", "r"):
  if "Best nsites:" in line:
    print line
    #line = line_nsites
    nsitesSplits = line.split(" ")
    nsitesFastep = int(nsitesSplits[6])
    print nsitesFastep

# Parsing seq file
seq = open("m1sea.pir","r")
seq = seq.read()
NumResInSeq = len(seq)
seq_first_char = seq[0].split()
seqfile = "m1sea.pir"

# Determining if S-SAD
atom = ['S', 'Se', 'Zn']
atom = atom[1]
if atom[0]:
  # Count nsites from seq
  if seq_first_char == ['M']: 
    nsitesSeq = seq.count('M') - 1
    print nsitesSeq
  else:
    nsitesSeq = seq.count('M')
    print nsitesSeq
elif atom[1]:
  # Count nsites from seq
  if seq_first_char == ['M']:
    nsitesSeq = seq.count('M') -1 + seq.count('C')
  else:
    nsitesSeq = seq.count('M') + seq.count('C')
else:
  # If atom is something else nsitesSeq=nsitesFastep
  nsitesSeq = nsitesFastep
    
# Compare seq sites and fastep sites
def asu(nsitesSeq, nsitesFastep):
  return round(nsitesSeq / nsitesFastep)
  
if atom[0]:
  asu_mol_Fastep = asu(nsitesSeq, nsitesFastep) # Calling the function
elif atom[1]:
  asu_mol_Fastep = asu(nsitesSeq, nsitesFastep)
else:
  asu_mol_Fastep = asu(nsitesSeq, nsitesFastep)
print ("Molecules per ASU (fast_ep):  "), asu_mol_Fastep

# Parsing MATTHEWS_COEF.xml
if(os.path.isfile('MATTHEWS_COEF.xml')):
  xmldoc_matthews = ET.parse("MATTHEWS_COEF.xml")
  root_matthews = xmldoc_matthews.getroot()
  for result_all in root_matthews.findall("result"):
     asu_mol_mat = float(result_all.get("nmol_in_asu"))
     print ("Molecules per ASU (Matthews): "), asu_mol_mat
     solvent = float(result_all.get("percent_solvent"))/100
     print ("Percent solvent: "), solvent

# Solvent and hand check
if solvent < 0.275:
  shelxe = open(os.path.join("0.25", "shelxe.sh"), "w")
  if "original" in fast_ep_string:
    shelxe_sh = ["shelxe sad sad_fa -s%.2f -a2 -m50 -b\n" % solvent,
           "f2mtz HKLIN sad.phs HKLOUT phs.mtz <<EOF\n",
           "CELL %s %s %s %s %s %s\n" % (a_cell, b_cell, 
           c_cell, alpha, beta, gamma),
           "SYMM %s\n" % SG,
           "LABOUT H K L F FOM PHI SIGF\n",
           "CTYP H H H F W P Q\n",
           "END\n" "EOF"]
    shelxe.writelines(shelxe_sh)
  else:
    shelxe_sh = ["shelxe sad sad_fa -s%.2f -a2 -m50 -b -i\n" % solvent,
           "f2mtz HKLIN sad.phs HKLOUT phs.mtz <<EOF\n",
           "CELL %s %s %s %s %s %s\n" % (a_cell, b_cell, 
           c_cell, alpha, beta, gamma),
           "SYMM %s\n" % SG,
           "LABOUT H K L F FOM PHI SIGF\n",
           "CTYP H H H F W P Q\n",
           "END\n" "EOF\n"]
    shelxe.writelines(shelxe_sh)
  shelxe.close()
  os.system("sh 0.25/shelxe.sh")
elif solvent >= 0.275 and solvent < 0.325:
  shelxe = open(os.path.join("0.30", "shelxe.sh"), "w")
  if "original" in fast_ep_string:
    shelxe_sh = ["shelxe sad sad_fa -s%.2f -a2 -m50 -b\n" % solvent,
           "f2mtz HKLIN sad.phs HKLOUT phs.mtz <<EOF\n",
           "CELL %s %s %s %s %s %s\n" % (a_cell, b_cell, 
           c_cell, alpha, beta, gamma),
           "SYMM %s\n" % SG,
           "LABOUT H K L F FOM PHI SIGF\n",
           "CTYP H H H F W P Q\n",
           "END\n" "EOF"]
    shelxe.writelines(shelxe_sh)
  else:
    shelxe_sh = ["shelxe sad sad_fa -s%.2f -a2 -m50 -b -i\n" % solvent,
           "f2mtz HKLIN sad.phs HKLOUT phs.mtz <<EOF\n",
           "CELL %s %s %s %s %s %s\n" % (a_cell, b_cell, 
           c_cell, alpha, beta, gamma),
           "SYMM %s\n" % SG,
           "LABOUT H K L F FOM PHI SIGF\n",
           "CTYP H H H F W P Q\n",
           "END\n" "EOF"]
    shelxe.writelines(shelxe_sh)
  shelxe.close()
  os.system("sh 0.30/shelxe.sh")
elif solvent >= 0.325 and solvent < 0.375:
  shelxe = open(os.path.join("0.35", "shelxe.sh"), "w")
  if "original" in fast_ep_string:
    shelxe_sh = ["shelxe sad sad_fa -s%.2f -a2 -m50 -b\n" % solvent,
           "f2mtz HKLIN sad.phs HKLOUT phs.mtz <<EOF\n",
           "CELL %s %s %s %s %s %s\n" % (a_cell, b_cell, 
           c_cell, alpha, beta, gamma),
           "SYMM %s\n" % SG,
           "LABOUT H K L F FOM PHI SIGF\n",
           "CTYP H H H F W P Q\n",
           "END\n" "EOF"]
    shelxe.writelines(shelxe_sh)
  else:
    shelxe_sh = ["shelxe sad sad_fa -s%.2f -a2 -m50 -b -i\n" % solvent,
           "f2mtz HKLIN sad.phs HKLOUT phs.mtz <<EOF\n",
           "CELL %s %s %s %s %s %s\n" % (a_cell, b_cell, 
           c_cell, alpha, beta, gamma),
           "SYMM %s\n" % SG,
           "LABOUT H K L F FOM PHI SIGF\n",
           "CTYP H H H F W P Q\n",
           "END\n" "EOF"]
    shelxe.writelines(shelxe_sh)
  shelxe.close()
  os.system("sh 0.35/shelxe.sh")
elif solvent >= 0.375 and solvent < 0.425:
  shelxe = open(os.path.join("0.40", "shelxe.sh"), "w")
  if "original" in fast_ep_string:
    shelxe_sh = ["shelxe sad sad_fa -s%.2f -a2 -m50 -b\n" % solvent,
           "f2mtz HKLIN sad.phs HKLOUT phs.mtz <<EOF\n",
           "CELL %s %s %s %s %s %s\n" % (a_cell, b_cell, 
           c_cell, alpha, beta, gamma),
           "SYMM %s\n" % SG,
           "LABOUT H K L F FOM PHI SIGF\n",
           "CTYP H H H F W P Q\n",
           "END\n" "EOF"]
    shelxe.writelines(shelxe_sh)
  else:
    shelxe_sh = ["shelxe sad sad_fa -s%.2f -a2 -m50 -b -i\n" % solvent,
           "f2mtz HKLIN sad.phs HKLOUT phs.mtz <<EOF\n",
           "CELL %s %s %s %s %s %s\n" % (a_cell, b_cell, 
           c_cell, alpha, beta, gamma),
           "SYMM %s\n" % SG,
           "LABOUT H K L F FOM PHI SIGF\n",
           "CTYP H H H F W P Q\n",
           "END\n" "EOF"]
    shelxe.writelines(shelxe_sh)
  shelxe.close()
  os.system("sh 0.40/shelxe.sh")
elif solvent >= 0.425 and solvent < 0.475:
  shelxe = open(os.path.join("0.45", "shelxe.sh"), "w")
  if "original" in fast_ep_string:
    shelxe_sh = ["shelxe sad sad_fa -s%.2f -a2 -m50 -b\n" % solvent,
           "f2mtz HKLIN sad.phs HKLOUT phs.mtz <<EOF\n",
           "CELL %s %s %s %s %s %s\n" % (a_cell, b_cell, 
           c_cell, alpha, beta, gamma),
           "SYMM %s\n" % SG,
           "LABOUT H K L F FOM PHI SIGF\n",
           "CTYP H H H F W P Q\n",
           "END\n" "EOF"]
    shelxe.writelines(shelxe_sh)
  else:
    shelxe_sh = ["shelxe sad sad_fa -s%.2f -a2 -m50 -b -i\n" % solvent,
           "f2mtz HKLIN sad.phs HKLOUT phs.mtz <<EOF\n",
           "CELL %s %s %s %s %s %s\n" % (a_cell, b_cell, 
           c_cell, alpha, beta, gamma),
           "SYMM %s\n" % SG,
           "LABOUT H K L F FOM PHI SIGF\n",
           "CTYP H H H F W P Q\n",
           "END\n" "EOF"]
    shelxe.writelines(shelxe_sh)
  shelxe.close()
  os.system("sh 0.45/shelxe.sh")
elif solvent >= 0.475 and solvent < 0.525:
  shelxe = open(os.path.join("0.50", "shelxe.sh"), "w")
  if "original" in fast_ep_string:
    shelxe_sh = ["shelxe sad sad_fa -s%.2f -a2 -m50 -b\n" % solvent,
           "f2mtz HKLIN sad.phs HKLOUT phs.mtz <<EOF\n",
           "CELL %s %s %s %s %s %s\n" % (a_cell, b_cell, 
           c_cell, alpha, beta, gamma),
           "SYMM %s\n" % SG,
           "LABOUT H K L F FOM PHI SIGF\n",
           "CTYP H H H F W P Q\n",
           "END\n" "EOF"]
    shelxe.writelines(shelxe_sh)
  else:
    shelxe_sh = ["shelxe sad sad_fa -s%.2f -a2 -m50 -b -i\n" % solvent,
           "f2mtz HKLIN sad.phs HKLOUT phs.mtz <<EOF\n",
           "CELL %s %s %s %s %s %s\n" % (a_cell, b_cell, 
           c_cell, alpha, beta, gamma),
           "SYMM %s\n" % SG,
           "LABOUT H K L F FOM PHI SIGF\n",
           "CTYP H H H F W P Q\n",
           "END\n" "EOF"]
    shelxe.writelines(shelxe_sh)
  shelxe.close()
  os.system("sh 0.50/shelxe.sh")
elif solvent >= 0.525 and solvent < 0.575:
  shelxe = open(os.path.join("0.55", "shelxe.sh"), "w")
  if "original" in fast_ep_string:
    shelxe_sh = ["shelxe sad sad_fa -s%.2f -a2 -m50 -b\n" % solvent,
           "f2mtz HKLIN sad.phs HKLOUT phs.mtz <<EOF\n",
           "CELL %s %s %s %s %s %s\n" % (a_cell, b_cell, 
           c_cell, alpha, beta, gamma),
           "SYMM %s\n" % SG,
           "LABOUT H K L F FOM PHI SIGF\n",
           "CTYP H H H F W P Q\n",
           "END\n" "EOF"]
    shelxe.writelines(shelxe_sh)
  else:
    shelxe_sh = ["shelxe sad sad_fa -s%.2f -a2 -m50 -b -i\n" % solvent,
           "f2mtz HKLIN sad.phs HKLOUT phs.mtz <<EOF\n",
           "CELL %s %s %s %s %s %s\n" % (a_cell, b_cell, 
           c_cell, alpha, beta, gamma),
           "SYMM %s\n" % SG,
           "LABOUT H K L F FOM PHI SIGF\n",
           "CTYP H H H F W P Q\n",
           "END\n" "EOF"]
    shelxe.writelines(shelxe_sh)
  shelxe.close()
  os.system("sh 0.55/shelxe.sh")
elif solvent >= 0.575 and solvent < 0.625:
  shelxe = open(os.path.join("0.60", "shelxe.sh"), "w")
  if "original" in fast_ep_string:
    shelxe_sh = ["shelxe sad sad_fa -s%.2f -a2 -m50 -b\n" % solvent,
           "f2mtz HKLIN sad.phs HKLOUT phs.mtz <<EOF\n",
           "CELL %s %s %s %s %s %s\n" % (a_cell, b_cell, 
           c_cell, alpha, beta, gamma),
           "SYMM %s\n" % SG,
           "LABOUT H K L F FOM PHI SIGF\n",
           "CTYP H H H F W P Q\n",
           "END\n" "EOF"]
    shelxe.writelines(shelxe_sh)
  else:
    shelxe_sh = ["shelxe sad sad_fa -s%.2f -a2 -m50 -b -i\n" % solvent,
           "f2mtz HKLIN sad.phs HKLOUT phs.mtz <<EOF\n",
           "CELL %s %s %s %s %s %s\n" % (a_cell, b_cell, 
           c_cell, alpha, beta, gamma),
           "SYMM %s\n" % SG,
           "LABOUT H K L F FOM PHI SIGF\n",
           "CTYP H H H F W P Q\n",
           "END\n" "EOF"]
    shelxe.writelines(shelxe_sh)
    os.system("sh 0.60shelxe.sh")
  shelxe.close()
  os.system("sh 0.60/shelxe.sh")
elif solvent >= 0.625 and solvent < 0.675:
  shelxe = open(os.path.join("0.65", "shelxe.sh"), "w")
  if "original" in fast_ep_string:
    shelxe_sh = ["shelxe sad sad_fa -s%.2f -a2 -m50 -b\n" % solvent,
           "f2mtz HKLIN sad.phs HKLOUT phs.mtz <<EOF\n",
           "CELL %s %s %s %s %s %s\n" % (a_cell, b_cell, 
           c_cell, alpha, beta, gamma),
           "SYMM %s\n" % SG,
           "LABOUT H K L F FOM PHI SIGF\n",
           "CTYP H H H F W P Q\n",
           "END\n" "EOF"]
    shelxe.writelines(shelxe_sh)
  else:
    shelxe_sh = ["shelxe sad sad_fa -s%.2f -a2 -m50 -b -i\n" % solvent,
           "f2mtz HKLIN sad.phs HKLOUT phs.mtz <<EOF\n",
           "CELL %s %s %s %s %s %s\n" % (a_cell, b_cell, 
           c_cell, alpha, beta, gamma),
           "SYMM %s\n" % SG,
           "LABOUT H K L F FOM PHI SIGF\n",
           "CTYP H H H F W P Q\n",
           "END\n" "EOF"]
    shelxe.writelines(shelxe_sh)
  shelxe.close()
  os.system("sh 0.65/shelxe.sh")
elif solvent >= 0.675 and solvent < 0.725:
  shelxe = open(os.path.join("0.70", "shelxe.sh"), "w")
  if "original" in fast_ep_string:
    shelxe_sh = ["shelxe sad sad_fa -s%.2f -a2 -m50 -b\n" % solvent,
           "f2mtz HKLIN sad.phs HKLOUT phs.mtz <<EOF\n",
           "CELL %s %s %s %s %s %s\n" % (a_cell, b_cell, 
           c_cell, alpha, beta, gamma),
           "SYMM %s\n" % SG,
           "LABOUT H K L F FOM PHI SIGF\n",
           "CTYP H H H F W P Q\n",
           "END\n" "EOF"]
    shelxe.writelines(shelxe_sh)
  else:
    shelxe_sh = ["shelxe sad sad_fa -s%.2f -a2 -m50 -b -i\n" % solvent,
           "f2mtz HKLIN sad.phs HKLOUT phs.mtz <<EOF\n",
           "CELL %s %s %s %s %s %s\n" % (a_cell, b_cell, 
           c_cell, alpha, beta, gamma),
           "SYMM %s\n" % SG,
           "LABOUT H K L F FOM PHI SIGF\n",
           "CTYP H H H F W P Q\n",
           "END\n" "EOF"]
    shelxe.writelines(shelxe_sh)
  shelxe.close()
  os.system("sh 0.70/shelxe.sh")
elif solvent >= 0.725:
  shelxe = open(os.path.join("0.75", "shelxe.sh"), "w")
  if "original" in fast_ep_string:
    shelxe_sh = ["shelxe sad sad_fa -s%.2f -a2 -m50 -b\n" % solvent,
           "f2mtz HKLIN sad.phs HKLOUT phs.mtz <<EOF\n",
           "CELL %s %s %s %s %s %s\n" % (a_cell, b_cell, 
           c_cell, alpha, beta, gamma),
           "SYMM %s\n" % SG,
           "LABOUT H K L F FOM PHI SIGF\n",
           "CTYP H H H F W P Q\n",
           "END\n" "EOF"]
    shelxe.writelines(shelxe_sh)
  else:
    shelxe_sh = ["shelxe sad sad_fa -s%.2f -a2 -m50 -b -i\n" % solvent,
           "f2mtz HKLIN sad.phs HKLOUT phs.mtz <<EOF\n",
           "CELL %s %s %s %s %s %s\n" % (a_cell, b_cell, 
           c_cell, alpha, beta, gamma),
           "SYMM %s\n" % SG,
           "LABOUT H K L F FOM PHI SIGF\n",
           "CTYP H H H F W P Q\n",
           "END\n" "EOF"]
    shelxe.writelines(shelxe_sh)
  shelxe.close()
  os.system("sh 0.75/shelxe.sh")
