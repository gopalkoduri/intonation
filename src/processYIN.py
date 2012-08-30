#!/usr/bin/env python
import numpy as np

def justify(yinpitch, justinpitch):
	#print yinpitch, justinpitch
	centDiff = 1200*np.log2(1.0*yinpitch/justinpitch)
	allowedThresh = 100 #cents
	if abs(centDiff) < allowedThresh:
		#print "returned", yinpitch
		return yinpitch
	elif np.mod(abs(centDiff), 1200) < allowedThresh:
		if centDiff > 0:
			foldFactor = 1.0/(int(centDiff/1200)+1)
		else:
			foldFactor = (int(abs(centDiff)/1200)+1)
		#print "returned", yinpitch*foldFactor
		return yinpitch*foldFactor
	else:
		#print "returned -1"
		return -1

def processYIN(yinfile, justinfile, yinjustinfile):
	#yin = np.loadtxt("../features/pitch-yin/2c0f7c85-32b1-4b05-84fe-eeaaa70c0fa6.wav.txt", delimiter=",", dtype="float")
	#justin = np.loadtxt("../features/pitch-new/2c0f7c85-32b1-4b05-84fe-eeaaa70c0fa6/2c0f7c85-32b1-4b05-84fe-eeaaa70c0fa6-Justin.txt", delimiter="\t", dtype="float")
	yin = np.loadtxt(yinfile, delimiter=",", dtype="float")
	justin = np.loadtxt(justinfile, delimiter="\t", dtype="float")
	#rounding the time values to 3 digits to allow comparison of frames.
	yin[:,0] = yin[:,0]
	step = 0.0029024956342978405
	for indyin in xrange(len(yin)):
		tyin = yin[indyin, 0]
		indjustin = int(tyin/step)+1
		try:
			if (justin[indjustin, 1] <= 0) or (yin[indyin, 1] <= 0):
				yin[indyin, 1] = -1
				continue
		except:
			print "index issue: ", indjustin, len(justin)
			continue
		correctedPitch = justify(yin[indyin, 1], justin[indjustin, 1])
		yin[indyin, 1] = correctedPitch
	np.savetxt(yinjustinfile, yin, fmt="%f", delimiter="\t")
	#return yin
