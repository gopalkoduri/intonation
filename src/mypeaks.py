#!/usr/bin/env python

import essentia.standard as es
import numpy as np
import intonationLib as iL

def nextJI(curJI):
	JICents = np.array([-1200, -1089, -997, -885, -814, -702, -591, -499, -387, -316, -204, -112, 0, 111, 203, 315, 386, 498, 609, 701, 813, 884, 996, 1088, 1200, 1311, 1403, 1515, 1586, 1698, 1809, 1901, 2013, 2084, 2196, 2288, 2400, 2511, 2603, 2715, 2786, 2898, 3009, 3101, 3213, 3284, 3396, 3488, 3600])
	curIndex = np.where(JICents == curJI)
	if curIndex[0][0]+1 < len(JICents):
		return JICents[curIndex[0][0]+1]
	else:
		raise IndexError, "Increase number of octaves in JICents!"


def nearestJI(ET):
	JICents = np.array([-1200, -1089, -997, -885, -814, -702, -591, -499, -387, -316, -204, -112, 0, 111, 203, 315, 386, 498, 609, 701, 813, 884, 996, 1088, 1200, 1311, 1403, 1515, 1586, 1698, 1809, 1901, 2013, 2084, 2196, 2288, 2400, 2511, 2603, 2715, 2786, 2898, 3009, 3101, 3213, 3284, 3396, 3488, 3600])
	if ET < JICents[0]-25 or ET > JICents[-1]+25:
		raise IndexError, "Increase number of octaves in JICents!"

	index = iL.findNearestIndex(JICents, ET)
	return JICents[index]

def peakdetect(n, binCenters, method="JI", window=100, peakAmpThresh=0.00005, valleyThresh=0.00003):
	"""
	This function expects smoothed histogram (i.e., n).

	method can be JI/ET/hybrid.
	JI and ET methods do not use generic peak detection algorithm. They use intervals and window
	to pick up the local maximum, and later filter out irrelevant peaks by using empirically found
	thresholds. Hybrid approach first finds peaks using generic peak detection algo, then filters
	the peaks heuristically.
	
	firstCenter and lastCenter can be given in multiples of 100. If the method is JI, they will be 
	automatically adjusted to nearest interval.

	window refers to the cent range used while picking up the maxima.

	The method returns:
	{"peaks":[[peak positions], [peak amplitudes]], "valleys": [[valley positions], [valley amplitudes]]}
	"""
	data = zip(binCenters, n)
	binCenters = np.array(binCenters)
	firstCenter = (min(binCenters)+1.5*window)/window*window
	lastCenter = (max(binCenters)-window)/window*window


	if method == "slope" or method == "hybrid":
		peaks = {}
		peakInfo = iL.findPeaks(n, binCenters, smoothingFactor=0, lookahead=20, delta=valleyThresh, averageHist=True)

		#find correspondences between peaks and valleys, and set valleys are left and right Indices
		#see the other method(s) for clarity!

		peakData = peakInfo["peaks"]
		valleyData = peakInfo["valleys"]

		#print len(peakData[0]), len(peakData[1])
		for i in xrange(len(peakData[0])):
			nearestIndex = iL.findNearestIndex(valleyData[0], peakData[0][i])
			if valleyData[0][nearestIndex] < peakData[0][i]:
				leftIndex = iL.findNearestIndex(binCenters, valleyData[0][nearestIndex])
				if (len(valleyData[0][nearestIndex+1:]) == 0):
					rightIndex = iL.findNearestIndex(binCenters, peakData[0][i]+window/2.0)
				else:
					offset = nearestIndex+1
					nearestIndex = offset+iL.findNearestIndex(valleyData[0][offset:], peakData[0][i])
					rightIndex = iL.findNearestIndex(binCenters, valleyData[0][nearestIndex])
			else:
				rightIndex = iL.findNearestIndex(binCenters, valleyData[0][nearestIndex])
				if (len(valleyData[0][:nearestIndex]) == 0):
					leftIndex = iL.findNearestIndex(binCenters, peakData[0][i]-window/2.0)
				else:
					nearestIndex = iL.findNearestIndex(valleyData[0][:nearestIndex], peakData[0][i])
					leftIndex = iL.findNearestIndex(binCenters, valleyData[0][nearestIndex])

			pos = iL.findNearestIndex(binCenters, peakData[0][i])
			print binCenters[pos], peakData[1][i], binCenters[leftIndex], binCenters[rightIndex]
			peaks[pos] = [peakData[1][i], leftIndex, rightIndex]

			if method == "hybrid": slopePeaks = peaks
	
	if method == "JI" or method == "ET" or method == "hybrid":
		peaks = {}
		#Obtain max value per interval
		if method == "JI" or method == "hybrid":
			firstCenter = nearestJI(firstCenter)
			lastCenter = nearestJI(lastCenter)

		interval = firstCenter
		prevInterval = firstCenter-window
		#NOTE: All *intervals are in cents. *indices are of binCenters/n
		while interval <= lastCenter:
			if method == "ET":
				leftIndex = iL.findNearestIndex(binCenters, interval-window/2)
				rightIndex = iL.findNearestIndex(binCenters, interval+window/2)
				interval += window
			elif method == "JI" or method == "hybrid":
				leftIndex = iL.findNearestIndex(binCenters, (interval+prevInterval)/2.0)
				prevInterval = interval
				interval = nextJI(interval)
				rightIndex = iL.findNearestIndex(binCenters, (interval+prevInterval)/2.0)
			peakPos = np.argmax(n[leftIndex:rightIndex])
			peakAmp = n[leftIndex+peakPos]
			peaks[leftIndex+peakPos] = [peakAmp, leftIndex, rightIndex]
			
			#print binCenters[leftIndex], binCenters[rightIndex], binCenters[leftIndex+peakPos], peakAmp
			#NOTE: All the indices (left/rightIndex, peakPos) are to be changed to represent respective cent 
			#value corresponding to the bin. Right now, they are indices of respective binCenters in the array.
	
	if method == "hybrid":
		#Mix peaks from slope method and JI method.
		p1 = slopePeaks.keys()
		p2 = peaks.keys()
		allPeaks = {} #overwriting peaks dict
		for p in p1:
			nearIndex = iL.findNearestIndex(p2, p)
			if abs(p-p2[nearIndex]) < 20: p2.pop(nearIndex)
	
		for p in p1: allPeaks[p] = slopePeaks[p]
		for p in p2: allPeaks[p] = peaks[p]
		peaks = allPeaks

	#Filter the peaks and retain eligible peaks, also get their valley points.

	# ----> peakAmpThresh <---- : remove the peaks which are below that

	for pos in peaks.keys():
		#pos is an index in binCenters/n. DOES NOT refer to a cent value.
		if peaks[pos][0] < peakAmpThresh:
			print "peakAmp: ", binCenters[pos]
			peaks.pop(pos)

	#Check if either left or right valley is deeper than ----> valleyThresh <----.
	valleys = {}
	for pos in peaks.keys():
		leftLobe = n[peaks[pos][1]:pos]
		rightLobe = n[pos:peaks[pos][2]]
		#Sanity check: Is it a genuine peak? Size of distributions on either side of the peak should be comparable.
		if len(rightLobe)!= 0 and (1.0*len(leftLobe)/len(rightLobe) < 0.15 or 1.0*len(leftLobe)/len(rightLobe) > 6.67):
			print "size: ", binCenters[pos]
			peaks.pop(pos)
			continue

		leftValleyPos = np.argmin(leftLobe)
		rightValleyPos = np.argmin(rightLobe)
		if (abs(leftLobe[leftValleyPos]-n[pos]) < valleyThresh and abs(rightLobe[rightValleyPos]-n[pos]) < valleyThresh):
			print "valley: ", binCenters[pos]
			peaks.pop(pos)
		else:
			valleys[peaks[pos][1]+leftValleyPos] = leftLobe[leftValleyPos]
			valleys[pos+rightValleyPos] = rightLobe[rightValleyPos]
	
	temp1 = np.array(peaks.values())
	print temp1
	temp1 = temp1[:, 0]

	return {'peaks':[binCenters[peaks.keys()], temp1], 'valleys':[binCenters[valleys.keys()], valleys.values()]}

