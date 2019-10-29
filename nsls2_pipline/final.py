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

seqfile = "m1sea.pir"
NumResInSeq = 100
solvent = 0.35
nsitesSeq = 2 
nsitesFastep = 10
nmol_in_asu = 1

eptree = ET.parse('pipelineoutput.xml')
nsitesFastep = int(eptree.findall('./CompInData/NSites')[0].text)
solvent   = float(eptree.findall('./CompInData/SolventContent')[0].text)

print ("Number of sites is: "), nsitesFastep
print ("Solvent content:    "), solvent

def asu(nsitesSeq, nsitesFastep):
  return round(nsitesSeq / nsitesFastep)

# Adding PHI and FOM from phs.mtz to sad.mtz
if (os.path.isfile("phs.mtz")):
	cad = open("cad.sh", "w")
	cad_sh = ["cad HKLIN1 sad.mtz HKLIN2 phs.mtz HKLOUT cad.mtz <<EOF\n",
			  "LABIN  FILE 1 E1=F   E2=SIGF E3=FreeF_flag\n",
			  "LABOUT FILE 1 E1=F   E2=SIGF E3=FreeF_flag\n",
			  "CTYP   FILE 1 E1=F   E2=Q    E3=I\n",
			  "LABIN  FILE 2 E1=FOM E2=PHI\n",
			  "LABOUT FILE 2 E1=FOM E2=PHI\n",
			  "CTYP   FILE 2 E1=W   E2=P\n" "END\n" "EOF"]
	cad.writelines(cad_sh)
	cad.close()
	os.system("sh cad.sh")
#else:
#	print "CAD failed, couldn't add PHI and FOM to sad.mtz"

# Writes the sh file to run dm
if (os.path.isfile("cad.mtz")):
	dm = open("dm.sh", "w")
	dm_sh = ["dm HKLIN cad.mtz HKLOUT cad_dm.mtz << endd\n",
				"SOLC %s\n" % solvent,
				"MODE SOLV HIST MULT\n",
				"NCYCLE 40\n",
				"SCHEME ALL\n",
				"COMBINE PERT FREE 2\n",
				"LABIN FP=F SIGFP=SIGF PHIO=PHI FOMO=FOM FREE=FreeF_flag\n",
				"LABOUT PHIDM=PHIDM FOMDM=FOMDM\n",
				"END\n"
				"endd"]
	dm.writelines(dm_sh)	
	dm.close()
	os.system("sh dm.sh")
#else:
#	print "Density modification failed"
	
# Run arpWarp
if (os.path.isfile("cad_dm.mtz")):
	arp = open("arp.sh", "w")
	#Matt_out = open("MATTHEWS_COEF.xml", "r")
	#Matt_r = Matt_out.read()
	#MattDic = xmltodict.parse(Matt_r)
	nmol_in_asu = asu(nsitesSeq, nsitesFastep)
        print ("Number mol in ASU:  "), nmol_in_asu
	#asu = nmol_in_asu.replace(" ", "")
	arp_sh = ["auto_tracing.sh datafile $PWD/cad_dm.mtz ",
				"residues %s cgr %s seqin $PWD/%s " % (NumResInSeq, nmol_in_asu, seqfile),
				"modelin $PWD/sad.pdb workdir $PWD fp FP ",
				"sigfp SIGFP phibest PHI fom FOM ",
				"freelabin FreeF_flag buildingcycles 2",
				"parfile 1, restrref 1 rrcyc 1 resol "20 2" "]
	arp.writelines(arp_sh)
	arp.close()
	arp = open("arp.sh", "r")
	print arp.read()
	os.system("sh arp.sh")
#else:
#	print "Model buidling failed"
	









