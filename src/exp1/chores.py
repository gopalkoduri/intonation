#import yaml
#
#tdir = "/media/CompMusic/users/gkoduri/Tonic/candidates_phaseII/"
#minTonic = 100
#maxTonic = 250
#annotations = yaml.load(file("../annotations.yaml"))
#
#for mbid in mbids:
#    if mbid not in annotations.keys():
#        annotations[mbid] = {}
#    if "tonic" not in annotations[mbid].keys():
#        try:
#            data = yaml.load(file(tdir+mbid+"_0-180.txt"))
#        except:
#            print mbid
#            continue
#        predictions = data.split()
#        tonic = -1
#        if float(predictions[0]) < minTonic:
#            tonic = float(predictions[0])*2
#        elif float(predictions[0]) > maxTonic:
#            tonic = float(predictions[0])/2
#        else:
#            tonic = float(predictions[0])
#        annotations[mbid]["tonic"] = {'value':tonic, 'annotator':'Justin\'s script'}
#        print mbid, " : ", tonic

#from essentia.standard import *
#from numpy import *
#from essentia import Pool
#
#fileIn = YamlInput(filename="../features/segments/f39e2ca6-93ad-4760-8f22-ed922f0dfd73/10/000.sig")
#pool = fileIn()
#
#mfccMeanIndices = [6, 9, 10, 11] #to be removed
#mfccVarIndices = [11, 12] #to be kept
#tristimulusMeanIndices = [3] #to be removed
#
#ditchedKeys = ["spectral_flux.skew", "spectral_flux.mean", "tristimulus.var", "pitch_confidence.var", "spectral_rolloff.mean", "zero_crossing_rate.mean"]
#
#mfccMeans = list(pool["mfcc.mean"])
#newMfccMeans = [i for j, i in enumerate(mfccMeans) if j not in mfccMeanIndices]
#mfccVars = list(pool["mfcc.var"])
#newMfccVars = [i for j, i in enumerate(mfccVars) if j in mfccVarIndices]
#tristimulusMeans = list(pool["tristimulus.mean"])
#newtristimulusMeans = [i for j, i in enumerate(tristimulusMeans) if j not in tristimulusMeanIndices]
#
#pool.remove("mfcc.mean")
#pool.remove("mfcc.var")
#pool.remove("tristimulus.mean")
#
#for i in newMfccMeans:
#	pool.add("mfcc.mean", i)
#
#for i in newMfccVars:
#	pool.add("mfcc.var", i)
#
#for i in newtristimulusMeans:
#	pool.add("tristimulus.mean", i)
#
#for key in ditchedKeys:
#    pool.remove(key)
#    
#output = YamlOutput(filename="testDeletion.sig")
#output(pool)

#import yaml
#from os import listdir

#data = yaml.load(file("vocal-raaga-artist-constrained.yaml"))
#mbids = []
#
#for artist in data.keys():
#	for raaga in data[artist]['recordings'].keys():
#		for rec in data[artist]['recordings'][raaga]:
#			mbids.append(rec['uuid'])

#bq = listdir("../audio/bad-quality")
#short = listdir("../audio/short")
#mq = listdir("../audio/medium-quality")
#
#bad_mbids = []
#temp = [i.split('.')[0] for i in bq]
#bad_mbids.extend(temp)
#temp = [i.split('.')[0] for i in mq]
#bad_mbids.extend(temp)
#temp = [i.split('.')[0] for i in short]
#bad_mbids.extend(temp)
#
#for mbid in mbids:
#	if mbid not in bad_mbids:
#		print mbid

#mbids = []
#mbids.extend(listdir('../audio/laptop'))
#mbids.extend(listdir('../audio/lab'))
#mbids.extend(listdir('../audio/done'))


#for i in xrange(len(n)):
	#if n[i]!=0:
		#for j in xrange(int(log10(n[i]))):
			#logCents.append(i)
			
			

#import numpy as np
#from scikits.learn import mixture

#clusters = []
##sizes = [500, 3000, 7000, 8000, 15000, 800, 2000, 700, 500]
##means = [0, 390, 693, 883, 1197, 1580, 1895, 2100, 2400]

#means = [0, 390, 1200]

#for i in xrange(len(means)):
	##c = means[i]+np.random.standard_normal(sizes[i])*30
	#c = means[i]+np.random.standard_normal(500)*30
	#c = c.reshape(c.size, 1)
	#clusters.append(c)

#for c in clusters:
	#bins = int(max(c)-min(c))/2
	#hist(c, bins, alpha=0.7)

#obs = np.concatenate(clusters[0:])
#g = mixture.GMM(n_states = len(means))
#g.fit(obs)

#means = g.means
#covars = g.covars
#weights = g.weights

#for i in xrange(len(means)):
	#data = np.random.normal(means[i], covars[i], int(weights[i]*10000))
	#bins = int(max(data)-min(data))/2
	#hist(data, bins, histtype="step", alpha=0.7)

#from os import listdir
#import yaml
#import numpy as np

#mbids = listdir("../features/pitch-new/")
#ann = yaml.load(file("../annotations.yaml"))

#candDir = "/media/CompMusic/users/gkoduri/Tonic/candidates_phaseII/"
#candFile = file("exp1-tonic-candidates.txt", "w")
#tonicDir = "../tonic-annotation/"

#for mbid in mbids:
	#if mbid in ann.keys():
		#if "tonic" in ann[mbid]:
			#if ann[mbid]["tonic"]["annotator"] == "Justin's script":
				#try:
					#data = np.loadtxt(tonicDir+mbid+"_0-180.txt")
				#except:
					#print "No tonic info for", mbid
					#continue
				#ann[mbid]["tonic"] = {'value':float(data[0]), 'annotator':'gopal'}
				##data = data.astype(str)
				##candFile.write(mbid+"_0-180.wav\t")
				##candFile.write(" ".join(data[0][:5]))
				##candFile.write("\n")
		#else:
			#print "---> No tonic for", mbid
	#else:
		#print "--> No mbid in annotations"
#print ann
#yaml.dump(ann, file("../annotations.yaml", "w"))


#Cumulative Histogram

import yaml
from numpy import *
from os import listdir
from os.path import isfile

ann = yaml.load(file("/home/gopal/workspace/annotations.yaml"))
pitchDir = "/home/gopal/workspace/features/pitch-yinjustin/"
allPitch = []

for f in audioFiles:
	if not isfile("../audio/"+f):
		print f, "is not a file!"
		continue
	mbid = f[:-4]
	pitchFile = pitchDir+mbid+"/"+mbid+"-YINJustin.txt"
	pitch = loadtxt(pitchFile)
	
	valid_pitch = [0]*len(pitch)
	timeStep = 0.0029024956342978405
	try:
		vocalSegments = ann[mbid]['vocal']
		for segment in vocalSegments['segments']:
			start = int(segment['start']/timeStep)
			end = int(segment['end']/timeStep) - 1
			if end > len(pitch):
				end = len(pitch)
			i = start
			while (i < end):
				if pitch[i] > 0:
					valid_pitch[i] = pitch[i]
				i += 1
	except:
		valid_pitch = pitch
		print "No annotations found, using the full track of", mbid
	
	valid_pitch = [i for i in valid_pitch if i > 0]
	tonic = ann[mbid]['tonic']['value']
	cents = [1200*log2(i/tonic) for i in valid_pitch]
	allPitch.extend(cents)
	
#Retain Good MBIDs

#import yaml
#import numpy as np

#data = yaml.load(file("vocal-raaga-artist-constrained.yaml"))
#mbids = np.loadtxt("66mbids.txt", dtype="str")
#raagaMBIDs = {}

#count = 0
#for artistid in data.keys():
	#for raaga in data[artistid]["recordings"].keys():
		#if raaga not in raagaMBIDs.keys():
			#raagaMBIDs[raaga] = []
		#for rec in data[artistid]["recordings"][raaga]:
			#if rec["uuid"] not in mbids:
				##data[artistid]["recordings"][raaga].remove(rec)
				#print rec["uuid"], "is not good."
			#else:
				#count = count+1
				#raagaMBIDs[raaga].append(rec["uuid"])

#print count
#yaml.dump(data, file("goodMBIDs-02Apr.yaml", "w"), default_flow_style=False)
#yaml.dump(raagaMBIDs, file("raagaMBIDs.yaml", "w"), default_flow_style=False)


# Per raaga Histogram

#import yaml
#from numpy import *
#from os import listdir
#from os.path import isfile

#raagaMBIDs = yaml.load(file("raagaMBIDs.yaml"))
#raaga = "kamas"
#mbids = raagaMBIDs[raaga]

#ann = yaml.load(file("../annotations.yaml"))
#pitchDir = "../features/pitch-new/"
#allPitch = []

#for mbid in mbids:
	#pitchFile = pitchDir+mbid+"/"+mbid+"-YINJustin.txt"
	#if not isfile(pitchFile):
		#print pitchFile, "is not a file!"
		#continue
	#pitch = loadtxt(pitchFile)
	#print pitchFile
	
	#valid_pitch = [0]*len(pitch)
	#timeStep = 0.0029024956342978405
	#try:
		#vocalSegments = ann[mbid]['vocal']
		#for segment in vocalSegments['segments']:
			#start = int(segment['start']/timeStep)
			#end = int(segment['end']/timeStep) - 1
			#if end > len(pitch):
				#end = len(pitch)
			#i = start
			#while (i < end):
				#if pitch[i] > 0:
					#valid_pitch[i] = pitch[i]
				#i += 1
	#except:
		#valid_pitch = pitch
		#print "No annotations found, using the full track of", mbid
	
	#valid_pitch = [i for i in valid_pitch if i > 0]
	#tonic = ann[mbid]['tonic']['value']
	#cents = [1200*log2(i/tonic) for i in valid_pitch]
	#allPitch.extend(cents)

#savetxt(raaga+"-allPitch.txt", allPitch)
#def computerParameters(paths):
#	[kn, kbinCenters, kcents] = iL.computeHist(paths)
#	refPeakInfo = iL.findPeaks(kn, kbinCenters, lookahead=20, delta=0.00005, averageHist=True)
#	allParams = {}
#
#	for i in paths:
#		[n, binCenters, cents] = iL.computeHist([i])
#		peakInfo = iL.findPeaks(n, binCenters, lookahead=10, delta=0.00002, averageHist=False, refPeaks=refPeakInfo["peaks"])
#		params = iL.characterizePeaks(peakInfo["peaks"][0], peakInfo["valleys"][0], cents, distribThresh=50)
#		mbid = i.split("/")[-1][:-4]
#		allParams[mbid] = params
#	
#	return allParams
#
#for mbid in tonicAbsentPhaseII:
#	data = loadtxt("tonic-annotation/"+mbid+"_0-180.txt", dtype="str", delimiter="\t")
#	if mbid not in annotations.keys():
#		annotations[mbid.encode("UTF-8")] = {"tonic":{"value":float(data[0].encode("UTF-8")), "annotator":"gopal"}}
#	else:
#		annotations[mbid]["tonic"] = {"value":float(data[0].encode("UTF-8")), "annotator":"gopal"}