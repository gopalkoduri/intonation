
#import intonationLib as iL
#import DAO as dao
#from os import chdir, listdir
#from os.path import exists, abspath
#from os import mkdir

#chdir('Z:/users/gkoduri/workspace/')
#
#
##raagas = iL.raagasToFile("ragas.txt")
##analyzed = []
#
##for raaga in raagas:
#	#if raaga in analyzed: continue
#	#recordings = iL.recordingsByRaagaAliases(raaga)
#	#if recordings:
#		#fwriter = file("raagas/hindustani/"+raaga+".txt", "w+")
#		#for rec in recordings.values():
#			#if rec[1] == None:
#				#fwriter.write("None\t"+rec[0]+"\n")
#			#else:
#				#fwriter.write(rec[1]+"\t"+rec[0]+"\n")
#			#analyzed.append(rec[0])
#		#fwriter.close()
#		#analyzed = list(set(analyzed))
#
#audioDir = "audio-hindustani/3/"
#pitchDirPrefix = "features/pitch-justin/"
##dirSuffix = "/alapanas/"
##dirSuffix = "/kritis/"
#dirSuffix = "/"
#for audiofile in listdir(abspath(audioDir)):
#	if audiofile.endswith(".wav"):
#		pitchdir = pitchDirPrefix+audiofile[:-4]+dirSuffix
#		if not exists(pitchdir):
#	##			if not exists(pitchDirPrefix+raaga):
#	##				mkdir(pitchDirPrefix+raaga)
#			mkdir(pitchdir)
#		pitchfile = audiofile[:-4]
#		iL.pitchJustin(audioDir+audiofile, pitchdir+pitchfile+"-Justin.txt")
#		#iL.pitchYIN(audiodir+audiofile, pitchdir+pitchfile+"-YIN.txt")

#fileindex['dheerasankarabharanam'] = iL.recordingsByRaagaAliases('sankarabharanam')
#fileindex = {}
#for raaga in raagas:
	#raaga = raaga[:-4]
	#fileindex[raaga] = iL.recordingsByRaagaAliases(raaga)
	
#for raaga in fileindex.keys():
	#for mbid in fileindex[raaga].keys():
		#print "copy", fileindex[raaga][mbid][1], "audio/"+raaga+"/"+mbid+".mp3"
		#if not os.path.exists("audio/"+raaga):
				#os.mkdir("audio/"+raaga)
		#if fileindex[raaga][mbid][1]:
			#shutil.copy(fileindex[raaga][mbid][1], "audio/"+raaga+"/"+mbid+".mp3")


#define _boxplot(swara)
#swaras = ["Sa", "R2/G1", "G3", "P", "D2/N1"]
#parameter = "skew2"
#khamasParams = {}
#kalyaniParams = {}
#for i in swaras:
	#khamasParams[i]=[]
	#kalyaniParams[i]=[]

#swaras = ["Sa_", "R1_", "R2/G1_", "G2/R3_", "G3_", "M1_", "M2_", "P_", "D1_", "D2/N1_", "D3/N2_", "N3_", "Sa", "R1", "R2/G1", "G2/R3", "G3", "M1", "M2", "P", "D1", "D2/N1", "D3/N2", "N3", "Sa^", "R1^", "R2/G1^", "G2/R3^", "G3^", "M1^", "M2^", "P^", "D1^", "D2/N1^", "D3/N2^", "N3^", "Sa^^", "R1^^", "R2/G1^^", "G2/R3^^", "G3^^", "M1^^", "M2^^", "P^^", "D1^^", "D2/N1^^", "D3/N2^^", "N3^^"]
#swaras = ["Sa", "R1", "R2/G1", "G2/R3", "G3", "M1", "M2", "P", "D1", "D2/N1", "D3/N2", "N3", "Sa^", "R1^", "R2/G1^", "G2/R3^", "G3^", "M1^", "M2^", "P^", "D1^", "D2/N1^", "D3/N2^", "N3^"]
#others = ["Sa_", "R1_", "R2/G1_", "G2/R3_", "G3_", "M1_", "M2_", "P_", "D1_", "D2/N1_", "D3/N2_", "N3_", "Sa^^", "R1^^", "R2/G1^^", "G2/R3^^", "G3^^", "M1^^", "M2^^", "P^^", "D1^^", "D2/N1^^", "D3/N2^^", "N3^^"]
#
#for raaga in allParameters.keys():
#	for rec in allParameters[raaga].keys():
#		for swara in swaras:
#			if swara not in allParameters[raaga][rec].keys():
#				allParameters[raaga][rec][swara] = {"position":0, "mean":0, "variance":0, "skew1":0, "kurtosis":0, "amplitude":0, "skew2":0}
#		for swara in allParameters[raaga][rec].keys():
#			if swara in others:
#				allParameters[raaga][rec].pop(swara)
#
#wekapath = "experiment-1/weka/"
#for raaga in allParameters.keys():
#	for rec in allParameters[raaga].keys():
		#print allParameters[raaga][rec]
		#yaml.dump(allParameters[raaga][rec], file(wekapath+raaga+"/"+rec+".sig", "w"), default_flow_style=False)

#performerRecs = {}

#for rec in recs:
	#info = db.getRecordingInfo("mbid", rec)
	#if info["artist"]["name"] in performerRecs.keys():
		#performerRecs[info["artist"]["name"]].append(rec)
	#else:
		#performerRecs[info["artist"]["name"]] = [rec]

#wekapath = "experiment-1/weka/performers/"
#for rec in paramsbyrec.keys():
#	info = db.getRecordingInfo("mbid", rec)
#	if exists(wekapath+info["artist"]["name"]):
#		yaml.dump(paramsbyrec[rec], file(wekapath+info["artist"]["name"]+"/"+rec+".sig", "w"), default_flow_style=False)

#import yaml
#import numpy as np
#
#allParameters = yaml.load(file("allParameters-25-50.yaml"))
#positions = {}
#allPositions = []
#for raaga in allParameters.keys():
#	positions[raaga] = {}
#	for rec in allParameters[raaga].keys():
#		offset = 0
#		if "Sa^" in allParameters[raaga][rec].keys():
#			offset = allParameters[raaga][rec]["Sa^"]-1200
#		for swara in allParameters[raaga][rec].keys():
#			if swara not in positions[raaga].keys():
#				positions[raaga][swara] = [allParameters[raaga][rec][swara]["mean"]-offset]
#			else:
#				positions[raaga][swara].append(allParameters[raaga][rec][swara]["mean"]-offset)
#			allPositions.append(allParameters[raaga][rec][swara]["mean"]-offset)
#
#from scipy.stats import mode, variation, skew
#allPositions = np.array(allPositions)
#equi = []
#means = []
#for i in xrange(24):
#	position = i*100
#	peak = allPositions[allPositions>position-50]
#	peak = peak[peak<position+50]
#	if len(peak) > 15:
#		print position, "\t", len(peak), "\t", np.mean(peak), "\t", mode(peak), "\t", variation(peak), "\t", skew(peak)
#		equi.append(position)
#		means.append(np.mean(peak))

#Average histogram plot in paper.
import intonationLib as iL
import yaml
from scipy.ndimage.filters import gaussian_filter
from pylab import *

raaga = yaml.load(file("../data/exp1/mukhari.yaml"))
kpaths = ["../../features/pitch-yinjustin/"+i+".txt" for i in raaga.keys()]
[kaln, kalbinCenters, kalcents] = iL.computeHist(kpaths[:8])
kalnS = gaussian_filter(kaln, 5)
refPeakInfo = iL.findPeaks(kalnS, kalbinCenters, smoothingFactor=5, lookahead=20, delta=0.00005, averageHist=True)

alln = []
allbinCenters = []
allnS = []

for i in kpaths[:8]:
	[n, binCenters, cents] = iL.computeHist([i])
	alln.append(n)
	allbinCenters.append(binCenters)
	allnS.append(gaussian_filter(n, 5))

_len = len(allnS)
for i in xrange(_len):
	plot(allbinCenters[i], allnS[i], lw=0.5)

plot(kalbinCenters, kalnS, lw=6, c="green")
xlabel("Pitch value (Cents)", fontsize=30)
ylabel("Normalized frequency", fontsize=30)
xticks(fontsize=22)
yticks(fontsize=22)
plot(refPeakInfo["peaks"][0], refPeakInfo["peaks"][1], 'rD', ms=10)
plot(refPeakInfo["valleys"][0], refPeakInfo["valleys"][1], 'yD', ms=10)
