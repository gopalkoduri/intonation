# -*- coding: utf-8 -*-
"""
If you don't see the font below or see question marks or boxes,
get the following unicode Telugu font and install.

Font: http://web.nickshanks.com/downloads/fonts/Pothana.ttf

@author: Gopala Krishna Koduri/గోపాల కృష్ణ కోడూరి
"""

from os import listdir
import essentia
from essentia.standard import *
from pylab import *
import numpy
import re
import sys

def sortNumerical(s):
	return int(s[:-4])

def toWeka(pathGiven, label, descriptor_flag=1):
	#assuming pathGiven will be like "so/on/and/so/forth/blah/blah2/raaga-name/"
	parts = pathGiven.strip("/").split("/")
	parentDir = "/".join(parts[:-1])
	arffFilename = "/".join([parentDir, parts[-1]+".arff"])
	wekafile=file(arffFilename, "w+")
	if descriptor_flag == 0:
		relation_name = "classification"
		wekafile.write("@RELATION classification\n\n")
		all_labels = ["ananda-bhairavi", "bhairavi", "hindolam", "kalyani", "khamas", "pantuvarali", "saveri", "sourashtram", "thodi"]
		#all_labels = ['Abhishek-Raghuram', 'S.-Sowmya', 'Aruna-Sairam', 'D.-K.-Pattammal', 'D.-K.-Jayaraman', 'T.-N.-Seshagopalan', 'Sanjay-Subrahmanyan', 'P.-Unnikrishnan', 'T.-M.-Krishna', 'Ranjani-Gayatri', 'G.-N.-Balasubramaniam', 'M.-D.-Ramanathan', 'Dandapani-Desikar']
		
	extension = ".sig"
	fnames = listdir(pathGiven)
	#fnames.sort(key=sortNumerical)
	for fname in fnames:
		if (fname.endswith(extension)):
			fname = pathGiven+"/"+fname
			print fname
			yamlinput = YamlInput(filename=fname)
			pool = yamlinput()
			descriptorList = pool.descriptorNames()
			#descriptorList.remove('metadata.version.essentia')
			
			#write descriptors details (once)
			if descriptor_flag == 0:
				attr_list = ""
				for i in descriptorList:
					if isinstance(pool[i], float):
						attr_list = attr_list + "@ATTRIBUTE "+i+" REAL\n"
					elif isinstance(pool[i], numpy.ndarray):
						for j in xrange(len(pool[i])):
							attr_list = attr_list + "@ATTRIBUTE "+i+str(j+1)+" REAL\n"
				attr_list = attr_list+"\n@ATTRIBUTE segment {"+", ".join(all_labels)+"}\n\n@DATA\n"
				wekafile.write(attr_list)
				descriptor_flag = 1
			
			#write the data points
			data_entry = ""
			for i in descriptorList:
				if isinstance(pool[i], float):
					data_entry = data_entry+str(pool[i])+", "
				elif isinstance(pool[i], numpy.ndarray):
					data_entry = data_entry+", ".join(str(x) for x in pool[i])+", "
			data_entry = data_entry+label
			wekafile.write(data_entry+"\n")

if __name__ == '__main__':
	#For building model
	#toWeka('weka/ananda-bhairavi/', 'ananda-bhairavi', descriptor_flag=0)
	#toWeka('weka/bhairavi/', 'bhairavi', descriptor_flag=1)
	#toWeka('weka/hindolam/', 'hindolam', descriptor_flag=0)
	#toWeka('weka/kalyani/', 'kalyani', descriptor_flag=0)
	toWeka('weka/khamas/', 'khamas', descriptor_flag=0)
	#toWeka('weka/pantuvarali/', 'pantuvarali', descriptor_flag=1)
	#toWeka('weka/saveri/', 'saveri', descriptor_flag=1)
	toWeka('weka/sourashtram/', 'sourashtram', descriptor_flag=1)
	#toWeka('weka/thodi/', 'thodi', descriptor_flag=1)
	
	#toWeka('weka/performers/Abhishek Raghuram/', 'Abhishek-Raghuram', descriptor_flag=0)
	#toWeka('weka/performers/Aruna Sairam/', 'Aruna Sairam', descriptor_flag=1)
	#toWeka('weka/performers/Dandapani Desikar/', 'Dandapani Desikar', descriptor_flag=1)
	#toWeka('weka/performers/D. K. Jayaraman/', 'D. K. Jayaraman', descriptor_flag=1)
	#toWeka('weka/performers/D. K. Pattammal/', 'D. K. Pattammal', descriptor_flag=1)
	#toWeka('weka/performers/G. N. Balasubramaniam/', 'G. N. Balasubramaniam', descriptor_flag=1)
	#toWeka('weka/performers/M. D. Ramanathan/', 'M. D. Ramanathan', descriptor_flag=1)
	#toWeka('weka/performers/P. Unnikrishnan/', 'P. Unnikrishnan', descriptor_flag=1)
	#toWeka('weka/performers/Ranjani-Gayatri/', 'Ranjani-Gayatri', descriptor_flag=1)
	#toWeka('weka/performers/Sanjay Subrahmanyan/', 'Sanjay-Subrahmanyan', descriptor_flag=0)
	#toWeka('weka/performers/S. Sowmya/', 'S.-Sowmya', descriptor_flag=1)
	#toWeka('weka/performers/T. M. Krishna/', 'T.-M.-Krishna', descriptor_flag=1)
	#toWeka('weka/performers/T. N. Seshagopalan/', 'T. N. Seshagopalan', descriptor_flag=1)
	#
	#For using the built model
	#pathGiven = sys.argv[1][]
	#pathsGiven = sys.argv[1:]
	#for pathGiven in pathsGiven:
	#	toWeka(pathGiven, "?")
	print "Done!"
