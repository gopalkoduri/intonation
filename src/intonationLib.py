#!/usr/bin/env python
# -*- coding:utf-8 -*-

import sys
sys.path.append("/home/gopal/workspace/stringDuplicates/src/")

import numpy as np
import pylab as p
import subprocess
import yaml
import wave

import essentia.standard as es
import DAO as dao
import peakDetect as pd
from stringDuplicates import stringDuplicates

from os import system, mkdir
from os.path import exists
from uuid import uuid1
from numpy import mod
from mpl_toolkits.axes_grid.anchored_artists import AnchoredText
from scipy.ndimage.filters import gaussian_filter, median_filter
from scipy.stats import skew, variation, kurtosis

def findNearestIndex(arr,value):
	"""For a given value, the function finds the nearest value 
	in the array and returns its index."""
	arr = np.array(arr)
	index=(np.abs(arr-value)).argmin()
	return index

def plotIntervals(n, binCenters):
	"""This function plots just intonation intervals given 
	the x, y of a plot"""
	justIntonation = [['Sa', 1.0], ['R1',16.0/15.0], ['R2/G1',9.0/8.0],\
	['G2/R3',6.0/5.0], ['G3',5.0/4.0], ['M1',4.0/3.0], ['M2',64.0/45.0],\
	['P',3.0/2.0], ['D1',8.0/5.0], ['D2/N1',5.0/3.0],\
	['D3/N2',16.0/9.0], ['N3',15.0/8.0]]
	JICents = []
	for interval in justIntonation:
		JICents.append([interval[0], round(1200.0*np.log2(interval[1]), 2)])
	
	#Intervals
	spacing = 0.02*max(n)
	for interval in JICents:
		p.axvline(x=interval[1], ls='-.', c='g', lw='1')
		p.text(interval[1]+3, spacing+4*max(n)/5, interval[0], fontsize='x-small')
		if interval[1]-2400 >= min(binCenters):
			p.axvline(x=interval[1]-2400, ls='--', c='r', lw='0.5')
			p.text(interval[1]-2400+3, spacing+4*max(n)/5, interval[0], fontsize='x-small')
		if interval[1]-1200 >= min(binCenters):
			p.axvline(x=interval[1]-1200, ls=':', c='b', lw='0.5')
			p.text(interval[1]-1200+3, spacing+4*max(n)/5, interval[0], fontsize='x-small')
		if interval[1]+1200 <= max(binCenters):
			p.axvline(x=interval[1]+1200, ls=':', c='b', lw='0.5')
			p.text(interval[1]+1200+3, spacing+4*max(n)/5, interval[0], fontsize='x-small')
		if interval[1]+2400 <= max(binCenters):
			p.axvline(x=interval[1]+2400, ls='-.', c='r', lw='0.5')
			p.text(interval[1]+2400+3, spacing+4*max(n)/5, interval[0], fontsize='x-small')
		spacing = -1*spacing;

def plotHist(n, binCenters, peakInfo=None, smoothingFactor=5, newFigure=True):
	"""This function plots histogram together with its smoothed
	version and peak information if provided. Just intonation 
	intervals are plotted for a reference."""
	if newFigure:
		p.figure()
	p.plot(binCenters, n, ls='-', c='#00CC00', lw='1')
	nS = gaussian_filter(n, smoothingFactor)
	p.plot(binCenters, nS, ls='-', c='b', lw='1.5')
	
	justIntonation = [['Sa', 1.0], ['R1',16.0/15.0], ['R2/G1',9.0/8.0],\
	['G2/R3',6.0/5.0], ['G3',5.0/4.0], ['M1',4.0/3.0], ['M2',64.0/45.0],\
	['P',3.0/2.0], ['D1',8.0/5.0], ['D2/N1',5.0/3.0],\
	['D3/N2',16.0/9.0], ['N3',15.0/8.0]]
	JICents = []
	for interval in justIntonation:
		JICents.append([interval[0], round(1200.0*np.log2(interval[1]), 2)])
	
	##Recordings Info
	#if len(filepaths) == 1:
		#db = dao.DAO()
		#recordingInfo = db.getRecordingInfo('mbid', mbid)
		#tags = ""
		#for tag in recordingInfo['tags']:
			#tags = tags+tag['tag'].capitalize()+" "+tag['category']+", "
		#tags = tags.strip(", ")
		#info = recordingInfo['title']+"\n"+\
		#recordingInfo['artist']['name']+"\n"+tags+\
		#"\nLength: "+strftime('%H:%M:%S', gmtime(recordingInfo['length']/1000.0))+\
		#"\nTonic: "+str(tonic)+" Hz"

		#ax = p.subplot(111)
		#at = AnchoredText(info,
					#prop=dict(size=10), frameon=True,
					#loc=2,
					#)
		#at.patch.set_boxstyle("round,pad=0.,rounding_size=0.1")
		#ax.add_artist(at)

	#Intervals
	spacing = 0.02*max(n)
	for interval in JICents:
		p.axvline(x=interval[1], ls='-.', c='g', lw='1')
		p.text(interval[1]+3, spacing+4*max(n)/5, interval[0], fontsize='x-small')
		if interval[1]-2400 >= min(binCenters):
			p.axvline(x=interval[1]-2400, ls='--', c='r', lw='0.5')
			p.text(interval[1]-2400+3, spacing+4*max(n)/5, interval[0], fontsize='x-small')
		if interval[1]-1200 >= min(binCenters):
			p.axvline(x=interval[1]-1200, ls=':', c='b', lw='0.5')
			p.text(interval[1]-1200+3, spacing+4*max(n)/5, interval[0], fontsize='x-small')
		if interval[1]+1200 <= max(binCenters):
			p.axvline(x=interval[1]+1200, ls=':', c='b', lw='0.5')
			p.text(interval[1]+1200+3, spacing+4*max(n)/5, interval[0], fontsize='x-small')
		if interval[1]+2400 <= max(binCenters):
			p.axvline(x=interval[1]+2400, ls='-.', c='r', lw='0.5')
			p.text(interval[1]+2400+3, spacing+4*max(n)/5, interval[0], fontsize='x-small')
		spacing = -1*spacing;

	p.title("Tonic-aligned pitch histogram without folding octaves")
	p.xlabel("Pitch value (Cents)")
	p.ylabel("Normalized frequency of occurence")
	if peakInfo:
		p.plot(peakInfo["peaks"][0], peakInfo["peaks"][1], 'rD', ms=5)
		p.plot(peakInfo["valleys"][0], peakInfo["valleys"][1], 'yD', ms=5)
	p.show()

def computeHist(filepaths, tonic=None, bins=None, annotationFile="/home/gopal/workspace/annotations.yaml", normed=True, folded=False, vocal=True, quantized=False, weight="duration", slopeThresh=1500, durationThresh=None): 
	"""Computes histogram of given array of filepaths. If the filepaths
	has more than one path, the returned histogram is an averaged one."""
	allCents = []
	mbid = ""
	for filepath in filepaths:
		print filepath
		mbid = filepath.split('/')[-1][:-4]
		if not tonic:
			try:
				manualAnnotations = yaml.load(file(annotationFile))
			except:
				print "The annotations file is missing, or the path is wrong. Please override the default argument value if needed. Quitting."
				return
			try:
				tonic = manualAnnotations[mbid]['tonic']['value']
			except:
				print mbid
				print "This file is neither annotated with tonic information nor you have passed it as an argument. Please try again. Quitting."
				return
		try:
			if quantized:
				data = np.array(noteSegment(pitchFile=filepath, slopeThresh=slopeThresh))
				i = 1
				_len = len(data)-1
				data[:,1] = median_filter(data[:,1], 7)
				count = 0
				while(i < _len):
					if data[i][1] == -1:
						i=i+1
						continue
					if (data[i][1]-data[i-1][1] != 0) and (data[i+1][1]-data[i][1] == 0):
						start = i
						while (i < _len and data[i+1][1]-data[i][1] == 0):
							i += 1
						if (data[i][0]-data[start][0])*1000 < durationThresh:
							count += 1
							i+=1
							data[start:i, 1] = np.zeros(i-start)-1
					else:
						i+=1
				data[:,1] = median_filter(data[:,1], 3)
				return data
			else:
				data = np.loadtxt(filepath, delimiter="\t", dtype="float")
		except:
			print filepath, ": this filepath is not correct. Please check. Quitting."
			return

		#Histogram
		dataLen = len(data)
		valid_pitch = [0]*dataLen
		#try:
			#taanamSegments = manualAnnotations[mbid]['taanam']
			#for segment in taanamSegments['segments']:
				#start = findNearestIndex(data[:,0], segment['start'])
				#end = findNearestIndex(data[:,0], segment['end'])
				#if end > dataLen:
					#end = dataLen
				#data[start:end, 1] = -1
		#except:
			#print "Taanam segments not found for", mbid
		if vocal:
			try:
				vocalSegments = manualAnnotations[mbid]['vocal']
				for segment in vocalSegments['segments']:
					start = findNearestIndex(data[:,0], segment['start'])
					end = findNearestIndex(data[:,0], segment['end'])
					if end > dataLen:
						end = dataLen
					for pitch_value in data[start:end, 1]:
						if pitch_value > 0:
							valid_pitch.append(pitch_value)
			except:
				valid_pitch = data[:,1]
				print "No annotated information of vocal segments found, using the full track for", mbid
		else:
			valid_pitch = data[:,1]

		valid_pitch = [i for i in valid_pitch if i > 0]
		cents = [1200*np.log2(i/tonic*1.0) for i in valid_pitch]
		if folded:
			cents = map(lambda x: int(x%1200), cents)
		tonic = None
		allCents.extend(cents)
	
	allCents = median_filter(allCents, 5)
	if not bins:
		bins = max(allCents)-min(allCents)
	if weight == "duration":
		n, binEdges = p.histogram(allCents, bins, normed=normed)
		binCenters = 0.5*(binEdges[1:]+binEdges[:-1])
		return [n, binCenters, allCents]
	else: # weight == "instance"
		if bins != 12:
			print "Quantization with bins != 12 is not handled. Quitting."
			exit()
		n = {}
		i = 1
		while(i < len(allCents)-1):
			if (allCents[i]-allCents[i-1] != 0) and (allCents[i+1]-allCents[i] == 0):
				flag = 0
				for val in xrange(allCents[i]-3, allCents[i]+3, 1):
				    if val in n.keys():
					n[val] += 1
					flag = 1
					break
				if flag == 0:
				    n[allCents[i]] = 1
			i+=1
		n = n.items()
		n.sort(key=lambda x:x[0])
		n = np.array(n)
		m = []
		for i in n[:,1]:
			m.append(1.0*i/sum(n[:,1]))
		return [np.array(m), n[:,0], allCents]

def characterizePeaks(peaks, valleys, cents, maxDistribThresh=50, minDistribThresh=25):
	"""Given peaks and cent values, it returns the characteristics
	of each peak."""
	justIntonation = [['Sa', 1.0], ['R1',16.0/15.0], ['R2/G1',9.0/8.0],\
	['G2/R3',6.0/5.0], ['G3',5.0/4.0], ['M1',4.0/3.0], ['M2',64.0/45.0],\
	['P',3.0/2.0], ['D1',8.0/5.0], ['D2/N1',5.0/3.0],\
	['D3/N2',16.0/9.0], ['N3',15.0/8.0]]
	swaras = ["Sa_", "R1_", "R2/G1_", "G2/R3_", "G3_", "M1_", "M2_", "P_", "D1_", "D2/N1_", "D3/N2_", "N3_", "Sa", "R1", "R2/G1", "G2/R3", "G3", "M1", "M2", "P", "D1", "D2/N1", "D3/N2", "N3", "Sa^", "R1^", "R2/G1^", "G2/R3^", "G3^", "M1^", "M2^", "P^", "D1^", "D2/N1^", "D3/N2^", "N3^", "Sa^^", "R1^^", "R2/G1^^", "G2/R3^^", "G3^^", "M1^^", "M2^^", "P^^", "D1^^", "D2/N1^^", "D3/N2^^", "N3^^"]
	JICents = []
	for interval in justIntonation:
		p = round(1200.0*np.log2(interval[1]), 2)
		JICents.append(p)
		JICents.append(p-1200)
		JICents.append(p+1200)
		JICents.append(p+2400)
		
	JICents.sort()
	cents = np.array(cents)
	parameters = {}
	for i in xrange(len(peaks[0])):
		p = peaks[0][i]
		#Set left and right bounds of the distribution.
		leftBound = p-maxDistribThresh
		rightBound = p+maxDistribThresh
		nearestValleyIndex = findNearestIndex(valleys[0], p)
		if p > valleys[0][nearestValleyIndex]:
			leftBound = valleys[0][nearestValleyIndex]
			if len(valleys[0][nearestValleyIndex+1:]) == 0:
				rightBound = p+maxDistribThresh
			else:
				offset = nearestValleyIndex+1
				nearestValleyIndex = findNearestIndex(valleys[0][nearestValleyIndex+1:], p)
				rightBound = valleys[0][offset+nearestValleyIndex]
		else:
			rightBound = valleys[0][nearestValleyIndex]
			if len(valleys[0][:nearestValleyIndex]) == 0:
				leftBound = p-maxDistribThresh
			else:
				nearestValleyIndex = findNearestIndex(valleys[0][:nearestValleyIndex], p)
				leftBound = valleys[0][nearestValleyIndex]
		
		if p-leftBound < minDistribThresh:
			leftBound = p-minDistribThresh
		if rightBound-p < minDistribThresh:
			rightBound = p+minDistribThresh
			
		newThresh = 0
		if p-leftBound < rightBound-p:
			newThresh = p-leftBound
		else:
			newThresh = rightBound-p
		if newThresh < maxDistribThresh:
			leftBound = p-newThresh
			rightBound = p+newThresh
		else:
			if p-leftBound > maxDistribThresh:
				leftBound = p-maxDistribThresh
			if rightBound-p > maxDistribThresh:
				rightBound = p+maxDistribThresh
		
		#extract the distribution and estimate the parameters
		distribution = cents[cents>leftBound]
		distribution = distribution[distribution<rightBound]
		#print p, "\t", len(distribution), "\t", leftBound, "\t", rightBound
		swaraIndex = findNearestIndex(JICents, p)
		swara = swaras[swaraIndex]
		_mean = float(np.mean(distribution))
		_variance = float(variation(distribution))
		_skew = float(skew(distribution))
		_kurtosis = float(kurtosis(distribution))
		pearsonSkew = float(3.0*(_mean-p)/np.sqrt(abs(_variance)))
		parameters[swara] = {"position":float(p), "mean":_mean, "variance":_variance, "skew1":_skew, "kurtosis":_kurtosis, "amplitude":float(peaks[1][i]), "skew2":pearsonSkew}
	return parameters

def findPeaks(n, binCenters, smoothingFactor=5, lookahead=15, delta=0.00002, averageHist=False, refPeaks=None, refThresh=25):
	#In the average histogram, we get the refPeaks and refValleys which are then used to clean the peaks we get from a single performance. For average histogram, lookahead should be ~ 20, delta = 0.000025 to pick up only the prominent peaks. For histogram of a performance, lookahead should be ~=15, delta = 0.00003 to pick up even the little peaks. 
	nS = gaussian_filter(n, smoothingFactor)
	_max, _min = pd.peakdetect(nS, binCenters, lookahead, delta)
	xPeaks = [p[0] for p in _max]
	yPeaks = [p[1] for p in _max]
	xValleys = [p[0] for p in _min]
	yValleys = [p[1] for p in _min]
	if averageHist:
		refPeaks = [xPeaks, yPeaks]
		refValleys = [xValleys, yValleys]
		return {"peaks":refPeaks, "valleys":refValleys}
	else:
		if not refPeaks:
			print "Reference peaks are not provided. Quitting."
			return
		xCleanPeaks = []
		yCleanPeaks = []
		#octave propagation of the reference peaks
		propThresh = 30
		tempPeaks = [i+1200 for i in refPeaks[0]]
		tempPeaks.extend([i-1200 for i in refPeaks[0]])
		extendedPeaks = []
		extendedPeaks.extend(refPeaks[0])
		for i in tempPeaks:
			#if a peak exists around, don't add this new one.
			nearestInd = findNearestIndex(refPeaks[0], i)
			diff = refPeaks[0][nearestInd] - i
			diff = np.mod(abs(diff), 1200)
			if diff > propThresh:
				extendedPeaks.append(i)
		#print extendedPeaks
		for peakLocationIndex in xrange(len(xPeaks)):
			extPeakLocationIndex = findNearestIndex(extendedPeaks, xPeaks[peakLocationIndex])
			diff = xPeaks[peakLocationIndex] - extendedPeaks[extPeakLocationIndex]
			diff = np.mod(abs(diff), 1200)
			#print xPeaks[peakLocationIndex], extendedPeaks[extPeakLocationIndex], diff
			if diff < refThresh:
				xCleanPeaks.append(xPeaks[peakLocationIndex])
				yCleanPeaks.append(yPeaks[peakLocationIndex])
		return {"peaks":[xCleanPeaks, yCleanPeaks], "valleys":[xValleys, yValleys]}

def computeParameters(paths):
	"""This is a wrapper function to bulk calculate parameters for a
	large number of files (run for one raaga at a time)"""
	print "Computing Average Histogram ..."
	[kn, kbinCenters, kcents] = computeHist(paths)
	refPeakInfo = findPeaks(kn, kbinCenters, lookahead=20, delta=0.00005, averageHist=True)
	allParams = {}
	
	print "Computing Parameters ..."
	for i in paths:
		[n, binCenters, cents] = computeHist([i])
		peakInfo = findPeaks(n, binCenters, lookahead=20, delta=0.000025, averageHist=False, refPeaks=refPeakInfo["peaks"])
		params = characterizePeaks(peakInfo["peaks"], peakInfo["valleys"], cents, maxDistribThresh=50, minDistribThresh=25)
		mbid = i.split("/")[-1][:-4]
		allParams[mbid] = params

	return allParams

def plotPitch(filepath, annotationFile="/home/gopal/workspace/annotations.yaml"):
	"""
	plotPitch(filepath):

	It plots the pitch contour based on the pitch information from the given file.
	The argument filepath refers to the file which has pitch information generated
	from Justin's code.
	"""

	data = file(filepath).readlines()
	time_stamps = []
	pitch_values = []
	mbid = filepath.split('/')[-1][:-4]
	for i in data:
		pitch_value = float(i.split("\t")[1].strip())
		time_stamps.append(float(i.split("\t")[0].strip()))
		if pitch_value > 0:
			pitch_values.append(pitch_value)
		else:
			pitch_values.append(0)
	#p.plot(time_stamps, pitch_values, ls="-", c='g', lw='0.5')
	#p.show();

	valid_pitch = [0]*len(pitch_values)
	timeStep = time_stamps[1] - time_stamps[0]
	try:
		manualAnnotations = yaml.load(file(annotationFile))
	except:
		print "I tried to read the annotation file to see if there are annotated vocal segments. But could not read the file. So I return the unfiltered pitch values"
		return [time_stamps, pitch_values]
	try:
		vocalSegments = manualAnnotations[mbid]['vocal']
		for segment in vocalSegments['segments']:
			for i in xrange(segment['start']/timeStep, segment['end']/timeStep-timeStep, timeStep):
				if pitch_values[i] > 0:
					valid_pitch[i] = pitch_values[i]
		print "I was able to read few annotated vocal segments, the pitch values in other than these segments are zeroed!"
		return [time_stamps, valid_pitch]
	except:
		print "I tried to read annotated vocal segments from the annotation file. But did not find the annotations in the file. So I return the unfiltered pitch values!"
		return [time_stamps, pitch_values]

def raagasToFile(filepath):
	"""
	raagasToFile(filepath)

	Writes all the raaga names from the database to the given file.
	"""
	db = dao.DAO()
	data = db.getTagsByCategory("raaga")
	ids = [data[1][i]['id'] for i in xrange(len(data[1]))]
	raagas = [data[1][i]['tag'] for i in xrange(len(data[1]))]
	#data = cursor.fetchall()
	filewriter = file(filepath, "w+")
	for i in xrange(len(raagas)):
		filewriter.write(str(ids[i])+"\t"+str(raagas[i])+"\n")
	filewriter.close()
	return raagas

def urlsByRaaga(raaga):
	"""
	urlsByRaaga(raaga)

	Returns musicbrainz urls of, recordings of a given raaga.
	"""
	db = dao.DAO()
	recordings = db.getRecordingsByTag("tag", raaga, category="raaga")
	if recordings[0]:
		uuids = [recordings[1][i]['uuid'] for i in xrange(len(recordings[1]))]
	else:
		print "No recordings found for", raaga, "raaga!!"
		return

	base_url="http://musicbrainz.org/recording/"
	urls = [base_url+mbid for mbid in uuids]
	return urls

def filepathByMBID(mbid, pathfile="/media/CompMusic/Carnatic/metadata/Carnatic.yaml"):
	"""
	filepathsByRaaga(raaga, pathfile="/media/CompMusic/Carnatic/metadata/Carnatic.yaml")

	Returns filepaths of recordings which are of given raaga.
	"""
	data = yaml.load(file(pathfile))
	if mbid in data.keys():
		return data[mbid]["path"]
	else:
		return None

def filepathsByRaaga(raaga):
	"""
	filepathsByRaaga(raaga)

	Returns filepaths of recordings which are of given raaga.
	"""
	db = dao.DAO()
	recordings = db.getRecordingsByTag("tag", raaga, category="raaga")
	fileindex = {}
	if recordings[0]:
		uuids = [recordings[1][i]['uuid'] for i in xrange(len(recordings[1]))]
	else:
		print "No recordings found for", raaga, "raaga!!"
		return

	for mbid in uuids:
		path = filepathByMBID(mbid)
		fileindex[mbid] = path

	return fileindex #fileindex[mbid]=path

def recordingsByRaagaAliases(raaga, sim_thresh=0.7):
	db = dao.DAO()
	tags = db.getTagsByCategory("raaga")
	raagas = [i["tag"] for i in tags[1]]
	aliases = stringDuplicates(raaga, raagas, sim_thresh)
	recordings = {}
	#recordings[mbid] = path

	for alias in aliases:
		data = db.getRecordingsByTag("tag", alias, category = "raaga")
		if data[0]:
			for i in xrange(len(data[1])):
				uuid = data[1][i]['uuid']
				path = filepathByMBID(uuid)
				recordings[uuid] = [alias, path]
	return recordings

def pitchJustin(audiofile, pitchfile, chunksize=300, overlap=0.29024956342978405, executable = "Z:/users/gkoduri/JustinPitchEvalBySynth/SalamonGomezMirex2011a/StreamEstimation.exe"):
	"""
	A function to extract pitch from longer files, using Justin's code. Input arguments
	are audio file's path, the chunksize, overlap (in seconds) and the path to executable.
	The overlap must be a multiple of 10 when converted to milliseconds.
	"""
	stepSize = 0.0029024956342978405 #NOTE: HARDCODED Values!!
	rem = mod(chunksize, stepSize)
	chunksize = chunksize-rem
	wfile = wave.open (audiofile, "r")
	length = (1.0 * wfile.getnframes ()) / wfile.getframerate ()
	if length > chunksize:
		start = 0
		finalFile = file(pitchfile, "w+")
		tempdir = "tmpPitch/"
		if not exists(tempdir):
			mkdir(tempdir)

		while start<length:
			base = tempdir+str(uuid1())
			newaudiofile = base+".wav"
			newpitchfile = base+".txt"
			end = start+chunksize
			if end > length:
				end = length
			loader = es.EasyLoader(filename=audiofile, startTime=start, endTime=end)
			audio = loader()
			writer = es.MonoWriter(filename = newaudiofile)
			writer(audio)
			command = executable+' "'+newaudiofile+'" "'+newpitchfile+'"'
			print command
			proc = subprocess.Popen(command).wait()
			nextStart = start+chunksize-overlap
			data = np.loadtxt(newpitchfile, dtype="float", delimiter="\t")
			data[:, 0] = data[:, 0]+start
			stepSize = int(overlap/chunksize)
			tobeDeleted = []
			for i in xrange(len(data[-1*stepSize:])):
				if data[i][0] >= nextStart:
					tobeDeleted.append(i)

			data = np.delete(data, tobeDeleted, axis=0)
			for i in xrange(len(data)):
					finalFile.write(str(data[i][0])+"\t"+str(data[i][1])+"\n")
			print start, tobeDeleted
			start=nextStart
		finalFile.close()
	else:
		command = executable+" "+audiofile+" "+pitchfile
		print command
		proc = subprocess.Popen(command).wait()

def pitchYIN(audiofile, pitchfile, framesize=2048, hopsize=1024):
	"""
	DO NOT use this!!
	"""
	pitch = es.PitchDetection(frameSize=framesize)
	spectrum = es.Spectrum(size=framesize)
	w = es.Windowing(size=framesize)
	loader = es.MonoLoader(filename=audiofile)
	audio = loader()

	pitch_values = []
	for frame in es.FrameGenerator(audio, framesize, hopsize):
		pitch_values.append(pitch(spectrum(w(frame))))

	reader = es.MetadataReader(filename=audiofile)
	metadata = reader()
	stepSize = 1.0*hopsize/metadata[9]

	fwriter = file(pitchfile, "w+")
	currentTime = 0
	for i in pitch_values:
		fwriter.write(str(round(currentTime, 6))+"\t"+str(round(i[0], 3))+"\n")
		currentTime = currentTime+stepSize
	fwriter.close()

def noteSegment(pitchFile, slopeThresh=1500, centsThresh=50, tonic=None, annotationFile="/home/gopal/workspace/annotations.yaml"):
	temp = pitchFile.split("/")
	f = temp[-1]
	checkPath = "/home/gopal/workspace/features/pitch-segmented-"+str(slopeThresh)+"/"+f
	if exists(checkPath):
		data = np.loadtxt(checkPath, dtype="float", delimiter="\t")
		return data
	
	data = np.loadtxt(pitchFile, dtype="float", delimiter="\t")
	eps = np.finfo(float).eps
	pitch = data[:,1]
	pitch = median_filter(pitch, 7)+eps
	time = data[:,0]
	justIntonation = [['Sa', 1.0], ['R1',16.0/15.0], ['R2/G1',9.0/8.0],\
	['G2/R3',6.0/5.0], ['G3',5.0/4.0], ['M1',4.0/3.0], ['M2',64.0/45.0],\
	['P',3.0/2.0], ['D1',8.0/5.0], ['D2/N1',5.0/3.0],\
	['D3/N2',16.0/9.0], ['N3',15.0/8.0]]
	JICents = []
	for interval in justIntonation:
		p = int(1200.0*np.log2(interval[1]))
		JICents.append(p)
		JICents.append(p-1200)
		JICents.append(p+1200)
		JICents.append(p+2400)
		
	JICents.sort()
	
	mbid = pitchFile.split('/')[-1][:-4]
	if not tonic:
		try:
			manualAnnotations = yaml.load(file(annotationFile))
		except:
			print "The annotations file is missing, or the path is wrong. Please override the default argument value if needed. Quitting."
			return
		try:
			tonic = manualAnnotations[mbid]['tonic']['value']
		except:
			print "This file is neither annotated with tonic information nor you have passed it as an argument. Please try again. Quitting."
			return
	
	cents = []
	for i in pitch:
		if i > 0:
			cents.append(1200*np.log2(1.0*i/tonic))
		else:
			cents.append(-10000)
		
	pitchSnapped = np.zeros(len(cents))
	for i in xrange(1, len(cents)-1):
		if cents[i] != -10000:
			ind = findNearestIndex(cents[i], JICents)
			centsDiff = abs(cents[i] - JICents[ind])
			slopeBack = abs((cents[i] - cents[i-1])/(time[i]-time[i-1]))
			slopeFront = abs((cents[i+1] - cents[i])/(time[i+1]-time[i]))
			if (slopeFront < slopeThresh or slopeBack < slopeThresh) and centsDiff <= centsThresh:
				pitchSnapped[i] = JICents[ind]
		else:
			pitchSnapped[i] = -10000
	
	for i in xrange(len(pitchSnapped)):
		if pitchSnapped[i] == -10000:
			pitchSnapped[i] = -1
		else:
			pitchSnapped[i] = tonic*pow(2, pitchSnapped[i]/1200.0)
	return zip(time, pitchSnapped)
