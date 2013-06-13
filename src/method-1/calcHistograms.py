#!/usr/bin/env python

import intonationLib as iL
import yaml
import pickle
from os import listdir
from os.path import exists

raagasPath = "/homedtic/gkoduri/Workspace/intonationLib/data/raagas/5/"
histogramsPath = "/homedtic/gkoduri/Workspace/intonationLib/data/method-1/"
pitchDir = "/homedtic/gkoduri/Workspace/features/pitch/"
overwrite = False

files = listdir(raagasPath)
for f in files:
	print f
	path = raagasPath+f
	mbids = yaml.load(file(path))
	mbids = mbids.keys()
	
	for mbid in mbids:
		print mbid
		if (exists(histogramsPath+"vocal-histograms/"+mbid+".pickle") or exists(histogramsPath+"non-vocal-histograms/"+mbid+".pickle")) and not overwrite:
			continue
		pitchFile = pitchDir+mbid+".txt"
		try:
			[n, binCenters, cents] = iL.computeHist(filepaths=[pitchFile])
			pickle.dump([n, binCenters, cents], file(histogramsPath+"histograms/"+mbid+".pickle", "w"))
		except:
			f = file("bad.txt", "a")
			f.write(mbid+"\n")
			f.close()
	print "-----------"
