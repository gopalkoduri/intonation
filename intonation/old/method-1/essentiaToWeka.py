# -*- coding: utf-8 -*-
"""
If you don't see the font below or see question marks or boxes,
get the following unicode Telugu font and install.

Font: http://web.nickshanks.com/downloads/fonts/Pothana.ttf

@author: Gopala Krishna Koduri/గోపాల కృష్ణ కోడూరి
"""

from os import listdir
import essentia.standard as es
import numpy as np
import sys

#write the descriptor details only once
descriptor_flag = 0

def sortNumerical(s):
	return int(s[:-4])

def toWeka(pathGiven, label):
	print pathGiven
	parts = pathGiven.rstrip("/").split("/")
	parentDir = "/".join(parts[:-1])
	arffFilename = "/".join([parentDir, parts[-1]+".arff"])
	wekafile=file(arffFilename, "w+")
	relation_name = "segmentation"
	global descriptor_flag
	if descriptor_flag == 0: wekafile.write("@RELATION segmentation\n\n")

	all_labels = ['anandabhairavi', 'asaveri', 'begada', 'bhairavi', 'dheerasankarabharanam', 'hamsadhvani', 'kaapi', 'kalyani', 'kamas', 'kambhoji', 'madhyamavathi', 'pantuvarali', 'saurashtram', 'thodi']
	
	extension = ".sig"
	fnames = listdir(pathGiven)
	#fnames.sort(key=sortNumerical)
	for fname in fnames:
		if (fname.endswith(extension)):
			fname = pathGiven+"/"+fname
			#print fname
			yamlinput = es.YamlInput(filename=fname)
			pool = yamlinput()
			descriptorList = pool.descriptorNames()
			
			#write descriptors details (once)
			if descriptor_flag == 0:
				attr_list = ""
				for i in descriptorList:
					if isinstance(pool[i], float):
						attr_list = attr_list + "@ATTRIBUTE "+i+" REAL\n"
					elif isinstance(pool[i], np.ndarray):
						shape = pool[i].shape
						_len = shape[0]
						if len(shape) > 1:
							_len = shape[0]*shape[1]
						for j in xrange(_len):
							attr_list = attr_list + "@ATTRIBUTE "+i+str(j+1)+" REAL\n"
				attr_list = attr_list+"\n@ATTRIBUTE segment {"+", ".join(all_labels)+"}\n\n@DATA\n"
				wekafile.write(attr_list)
				descriptor_flag = 1
			
			#write the data points
			data_entry = ""
			for i in descriptorList:
				if isinstance(pool[i], float):
					data_entry = data_entry+str(pool[i])+", "
				elif isinstance(pool[i], np.ndarray):
					data = pool[i]
					shape = data.shape
					if len(shape) > 1:
						data = data.reshape(shape[0]*shape[1])
					data_entry = data_entry+", ".join(str(x) for x in data)+", "
			data_entry = data_entry+label
			wekafile.write(data_entry+"\n")

if __name__ == '__main__':
	#For building model
	toWeka("/homedtic/gkoduri/Workspace/intonationLib/data/method-2/experiments/data/abhogi", "abhogi")
	toWeka("/homedtic/gkoduri/Workspace/intonationLib/data/method-2/experiments/data/anandabhairavi", "anandabhairavi")
	toWeka("/homedtic/gkoduri/Workspace/intonationLib/data/method-2/experiments/data/asaveri", "asaveri")
	toWeka("/homedtic/gkoduri/Workspace/intonationLib/data/method-2/experiments/data/atana", "atana")
	toWeka("/homedtic/gkoduri/Workspace/intonationLib/data/method-2/experiments/data/begada", "begada")
	toWeka("/homedtic/gkoduri/Workspace/intonationLib/data/method-2/experiments/data/behag", "behag")
	toWeka("/homedtic/gkoduri/Workspace/intonationLib/data/method-2/experiments/data/bhairavi", "bhairavi")
	toWeka("/homedtic/gkoduri/Workspace/intonationLib/data/method-2/experiments/data/bilahari", "bilahari")
	toWeka("/homedtic/gkoduri/Workspace/intonationLib/data/method-2/experiments/data/devagandhari", "devagandhari")
	toWeka("/homedtic/gkoduri/Workspace/intonationLib/data/method-2/experiments/data/devamanohari", "devamanohari")
	toWeka("/homedtic/gkoduri/Workspace/intonationLib/data/method-2/experiments/data/dhanashri", "dhanashri")
	toWeka("/homedtic/gkoduri/Workspace/intonationLib/data/method-2/experiments/data/dhanyasi", "dhanyasi")
	toWeka("/homedtic/gkoduri/Workspace/intonationLib/data/method-2/experiments/data/dheerasankarabharanam", "dheerasankarabharanam")
	toWeka("/homedtic/gkoduri/Workspace/intonationLib/data/method-2/experiments/data/hamsadhvani", "hamsadhvani")
	toWeka("/homedtic/gkoduri/Workspace/intonationLib/data/method-2/experiments/data/hari_kambhoji", "hari_kambhoji")
	toWeka("/homedtic/gkoduri/Workspace/intonationLib/data/method-2/experiments/data/hindolam", "hindolam")
	toWeka("/homedtic/gkoduri/Workspace/intonationLib/data/method-2/experiments/data/jaunpuri", "jaunpuri")
	toWeka("/homedtic/gkoduri/Workspace/intonationLib/data/method-2/experiments/data/kaapi", "kaapi")
	toWeka("/homedtic/gkoduri/Workspace/intonationLib/data/method-2/experiments/data/kalyani", "kalyani")
	toWeka("/homedtic/gkoduri/Workspace/intonationLib/data/method-2/experiments/data/kamas", "kamas")
	toWeka("/homedtic/gkoduri/Workspace/intonationLib/data/method-2/experiments/data/kambhoji", "kambhoji")
	toWeka("/homedtic/gkoduri/Workspace/intonationLib/data/method-2/experiments/data/karaharapriya", "karaharapriya")
	toWeka("/homedtic/gkoduri/Workspace/intonationLib/data/method-2/experiments/data/kedaragaula", "kedaragaula")
	toWeka("/homedtic/gkoduri/Workspace/intonationLib/data/method-2/experiments/data/madhyamavathi", "madhyamavathi")
	toWeka("/homedtic/gkoduri/Workspace/intonationLib/data/method-2/experiments/data/manji", "manji")
	toWeka("/homedtic/gkoduri/Workspace/intonationLib/data/method-2/experiments/data/mohanam", "mohanam")
	toWeka("/homedtic/gkoduri/Workspace/intonationLib/data/method-2/experiments/data/mukhari", "mukhari")
	toWeka("/homedtic/gkoduri/Workspace/intonationLib/data/method-2/experiments/data/nagasvarali", "nagasvarali")
	toWeka("/homedtic/gkoduri/Workspace/intonationLib/data/method-2/experiments/data/natakurinji", "natakurinji")
	toWeka("/homedtic/gkoduri/Workspace/intonationLib/data/method-2/experiments/data/pantuvarali", "pantuvarali")
	toWeka("/homedtic/gkoduri/Workspace/intonationLib/data/method-2/experiments/data/poorikalyani", "poorikalyani")
	toWeka("/homedtic/gkoduri/Workspace/intonationLib/data/method-2/experiments/data/ranjani", "ranjani")
	toWeka("/homedtic/gkoduri/Workspace/intonationLib/data/method-2/experiments/data/reethi_gowlai", "reethi_gowlai")
	toWeka("/homedtic/gkoduri/Workspace/intonationLib/data/method-2/experiments/data/sahana", "sahana")
	toWeka("/homedtic/gkoduri/Workspace/intonationLib/data/method-2/experiments/data/saurashtram", "saurashtram")
	toWeka("/homedtic/gkoduri/Workspace/intonationLib/data/method-2/experiments/data/senchurutti", "senchurutti")
	toWeka("/homedtic/gkoduri/Workspace/intonationLib/data/method-2/experiments/data/shanmukhapriya", "shanmukhapriya")
	toWeka("/homedtic/gkoduri/Workspace/intonationLib/data/method-2/experiments/data/shree_ranjani", "shree_ranjani")
	toWeka("/homedtic/gkoduri/Workspace/intonationLib/data/method-2/experiments/data/shuddha_saveri", "shuddha_saveri")
	toWeka("/homedtic/gkoduri/Workspace/intonationLib/data/method-2/experiments/data/sindhu_bhairavi", "sindhu_bhairavi")
	toWeka("/homedtic/gkoduri/Workspace/intonationLib/data/method-2/experiments/data/surati", "surati")
	toWeka("/homedtic/gkoduri/Workspace/intonationLib/data/method-2/experiments/data/thodi", "thodi")
	toWeka("/homedtic/gkoduri/Workspace/intonationLib/data/method-2/experiments/data/vachaspati", "vachaspati")
	toWeka("/homedtic/gkoduri/Workspace/intonationLib/data/method-2/experiments/data/vasanta", "vasanta")
	toWeka("/homedtic/gkoduri/Workspace/intonationLib/data/method-2/experiments/data/yadukula_kambhoji", "yadukula_kambhoji")

	#For using the built model
	#pathGiven = sys.argv[1]
	#pathsGiven = sys.argv[1:]
	#for pathGiven in pathsGiven:
	#	toWeka(pathGiven, "?")
	print "Done!"
