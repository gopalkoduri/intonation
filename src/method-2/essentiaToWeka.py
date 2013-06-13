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
	toWeka("/homedtic/gkoduri/Workspace/intonationLib/data/method-2/experiments/data-positionAmplitude/abhogi", "abhogi")
	toWeka("/homedtic/gkoduri/Workspace/intonationLib/data/method-2/experiments/data-positionAmplitude/anandabhairavi", "anandabhairavi")
	toWeka("/homedtic/gkoduri/Workspace/intonationLib/data/method-2/experiments/data-positionAmplitude/asaveri", "asaveri")
	toWeka("/homedtic/gkoduri/Workspace/intonationLib/data/method-2/experiments/data-positionAmplitude/atana", "atana")
	toWeka("/homedtic/gkoduri/Workspace/intonationLib/data/method-2/experiments/data-positionAmplitude/begada", "begada")
	toWeka("/homedtic/gkoduri/Workspace/intonationLib/data/method-2/experiments/data-positionAmplitude/behag", "behag")
	toWeka("/homedtic/gkoduri/Workspace/intonationLib/data/method-2/experiments/data-positionAmplitude/bhairavi", "bhairavi")
	toWeka("/homedtic/gkoduri/Workspace/intonationLib/data/method-2/experiments/data-positionAmplitude/bilahari", "bilahari")
	toWeka("/homedtic/gkoduri/Workspace/intonationLib/data/method-2/experiments/data-positionAmplitude/devagandhari", "devagandhari")
	toWeka("/homedtic/gkoduri/Workspace/intonationLib/data/method-2/experiments/data-positionAmplitude/devamanohari", "devamanohari")
	toWeka("/homedtic/gkoduri/Workspace/intonationLib/data/method-2/experiments/data-positionAmplitude/dhanashri", "dhanashri")
	toWeka("/homedtic/gkoduri/Workspace/intonationLib/data/method-2/experiments/data-positionAmplitude/dhanyasi", "dhanyasi")
	toWeka("/homedtic/gkoduri/Workspace/intonationLib/data/method-2/experiments/data-positionAmplitude/dheerasankarabharanam", "dheerasankarabharanam")
	toWeka("/homedtic/gkoduri/Workspace/intonationLib/data/method-2/experiments/data-positionAmplitude/hamsadhvani", "hamsadhvani")
	toWeka("/homedtic/gkoduri/Workspace/intonationLib/data/method-2/experiments/data-positionAmplitude/hari_kambhoji", "hari_kambhoji")
	toWeka("/homedtic/gkoduri/Workspace/intonationLib/data/method-2/experiments/data-positionAmplitude/hindolam", "hindolam")
	toWeka("/homedtic/gkoduri/Workspace/intonationLib/data/method-2/experiments/data-positionAmplitude/jaunpuri", "jaunpuri")
	toWeka("/homedtic/gkoduri/Workspace/intonationLib/data/method-2/experiments/data-positionAmplitude/kaapi", "kaapi")
	toWeka("/homedtic/gkoduri/Workspace/intonationLib/data/method-2/experiments/data-positionAmplitude/kalyani", "kalyani")
	toWeka("/homedtic/gkoduri/Workspace/intonationLib/data/method-2/experiments/data-positionAmplitude/kamas", "kamas")
	toWeka("/homedtic/gkoduri/Workspace/intonationLib/data/method-2/experiments/data-positionAmplitude/kambhoji", "kambhoji")
	toWeka("/homedtic/gkoduri/Workspace/intonationLib/data/method-2/experiments/data-positionAmplitude/karaharapriya", "karaharapriya")
	toWeka("/homedtic/gkoduri/Workspace/intonationLib/data/method-2/experiments/data-positionAmplitude/kedaragaula", "kedaragaula")
	toWeka("/homedtic/gkoduri/Workspace/intonationLib/data/method-2/experiments/data-positionAmplitude/madhyamavathi", "madhyamavathi")
	toWeka("/homedtic/gkoduri/Workspace/intonationLib/data/method-2/experiments/data-positionAmplitude/manji", "manji")
	toWeka("/homedtic/gkoduri/Workspace/intonationLib/data/method-2/experiments/data-positionAmplitude/mohanam", "mohanam")
	toWeka("/homedtic/gkoduri/Workspace/intonationLib/data/method-2/experiments/data-positionAmplitude/mukhari", "mukhari")
	toWeka("/homedtic/gkoduri/Workspace/intonationLib/data/method-2/experiments/data-positionAmplitude/nagasvarali", "nagasvarali")
	toWeka("/homedtic/gkoduri/Workspace/intonationLib/data/method-2/experiments/data-positionAmplitude/natakurinji", "natakurinji")
	toWeka("/homedtic/gkoduri/Workspace/intonationLib/data/method-2/experiments/data-positionAmplitude/pantuvarali", "pantuvarali")
	toWeka("/homedtic/gkoduri/Workspace/intonationLib/data/method-2/experiments/data-positionAmplitude/poorikalyani", "poorikalyani")
	toWeka("/homedtic/gkoduri/Workspace/intonationLib/data/method-2/experiments/data-positionAmplitude/ranjani", "ranjani")
	toWeka("/homedtic/gkoduri/Workspace/intonationLib/data/method-2/experiments/data-positionAmplitude/reethi_gowlai", "reethi_gowlai")
	toWeka("/homedtic/gkoduri/Workspace/intonationLib/data/method-2/experiments/data-positionAmplitude/sahana", "sahana")
	toWeka("/homedtic/gkoduri/Workspace/intonationLib/data/method-2/experiments/data-positionAmplitude/saurashtram", "saurashtram")
	toWeka("/homedtic/gkoduri/Workspace/intonationLib/data/method-2/experiments/data-positionAmplitude/senchurutti", "senchurutti")
	toWeka("/homedtic/gkoduri/Workspace/intonationLib/data/method-2/experiments/data-positionAmplitude/shanmukhapriya", "shanmukhapriya")
	toWeka("/homedtic/gkoduri/Workspace/intonationLib/data/method-2/experiments/data-positionAmplitude/shree_ranjani", "shree_ranjani")
	toWeka("/homedtic/gkoduri/Workspace/intonationLib/data/method-2/experiments/data-positionAmplitude/shuddha_saveri", "shuddha_saveri")
	toWeka("/homedtic/gkoduri/Workspace/intonationLib/data/method-2/experiments/data-positionAmplitude/sindhu_bhairavi", "sindhu_bhairavi")
	toWeka("/homedtic/gkoduri/Workspace/intonationLib/data/method-2/experiments/data-positionAmplitude/surati", "surati")
	toWeka("/homedtic/gkoduri/Workspace/intonationLib/data/method-2/experiments/data-positionAmplitude/thodi", "thodi")
	toWeka("/homedtic/gkoduri/Workspace/intonationLib/data/method-2/experiments/data-positionAmplitude/vachaspati", "vachaspati")
	toWeka("/homedtic/gkoduri/Workspace/intonationLib/data/method-2/experiments/data-positionAmplitude/vasanta", "vasanta")
	toWeka("/homedtic/gkoduri/Workspace/intonationLib/data/method-2/experiments/data-positionAmplitude/yadukula_kambhoji", "yadukula_kambhoji")

	#For using the built model
	#pathGiven = sys.argv[1]
	#pathsGiven = sys.argv[1:]
	#for pathGiven in pathsGiven:
	#	toWeka(pathGiven, "?")
	print "Done!"
