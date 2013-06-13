#!/usr/bin/env python
# -*- coding:utf-8 -*-
from __future__ import division

import sys
sys.path.append("/media/CompMusic/audio/users/gkoduri/Workspace/stringDuplicates/src/")

import numpy as np
import pylab as p
import subprocess
import yaml
import wave
import pickle

import essentia.standard as es
import DAO as dao
import peakDetection as peaks
from stringDuplicates import stringDuplicates

from os import system, mkdir
from os.path import exists
from copy import deepcopy
from uuid import uuid1
from numpy import mod
from mpl_toolkits.axes_grid.anchored_artists import AnchoredText
from scipy.ndimage.filters import gaussian_filter, median_filter
from scipy.stats import skew, variation, kurtosis


#Handy general purpose functions

def findNearestIndex(arr,value):
    """For a given value, the function finds the nearest value 
    in the array and returns its index."""
    arr = np.array(arr)
    index=(np.abs(arr-value)).argmin()
    return index

def rangeFilter(arr, _min, _max, scale="cents"):
    def _filter(x):
        if x >= _min and x <= _max:
            return x
        else:
            if scale == "cents":
                return -10000
            else:
                return -1
    res = [_filter(i) for i in arr]
    return res


#Functions for computation
def discreteFilter(mbid, srcdata=None, slopeThresh=1500, centsThresh=50, returnScale="cents", tonic=None, annotationFile="/media/CompMusic/audio/users/gkoduri/Workspace/annotations.yaml"):
    """If 'srcdata' is given, discreteFilter works on that. Otherwise, it works with mbid. Data, if given should be in
    cent scale!"""

    #checkPath = "/media/CompMusic/audio/users/gkoduri/Workspace/features/pitch-segmented-"+str(slopeThresh)+"/"+mbid+".txt"
    #if exists(checkPath):
    #    data = np.loadtxt(checkPath, dtype="float", delimiter="\t")
    
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
    
    #eps = np.finfo(float).eps
    #pitch = median_filter(pitch, 7)+eps
    if srcdata == None:
        data = pitch(mbid, returnScale="cents")
    else:
        data = deepcopy(srcdata)
    data[:, 1] =  median_filter(data[:, 1], 7)
    dataLen = len(data)
    pitchSnapped = np.zeros(dataLen)-1
    for i in xrange(1, dataLen-1):
        if data[i, 1] != -10000:
            ind = findNearestIndex(data[i, 1], JICents)
            centsDiff = abs(data[i, 1] - JICents[ind])
            slopeBack = abs((data[i, 1] - data[i-1, 1])/(data[i, 0]-data[i-1, 0]))
            slopeFront = abs((data[i+1, 1] - data[i, 1])/(data[i+1, 0]-data[i, 0]))
            if (slopeFront < slopeThresh or slopeBack < slopeThresh) and centsDiff <= centsThresh:
                if returnScale == "cents":
                    pitchSnapped[i] = JICents[ind]
                else:
                    pitchSnapped[i] = tonic*pow(2, JICents[ind]/1200.0)
        else:
            if returnScale == "cents":
                pitchSnapped[i] = -10000
            else:
                pitchSnapped[i] = -1
            
    
    data[:, 1] = pitchSnapped
    return data

def durationFilter(srcdata, durationThresh):
    """
    Gets rid of chunks of quantized pitches which last for less than
    a given threshold amount of time (in milliseconds). Input data must be
    a segmented/snapped pitch contour.
    """
    data = deepcopy(srcdata)
    i = 1
    _len = len(data)-1

    while(i < _len):
        if data[i][1] == -1:
            i=i+1
            continue
        if (data[i][1]-data[i-1][1] != 0) and (data[i+1][1]-data[i][1] == 0):
            start = i
            while (i < _len and data[i+1][1]-data[i][1] == 0):
                i += 1
            if (data[i][0]-data[start][0])*1000 < durationThresh:
                data[start:i+1, 1] = np.zeros(i+1-start)-1
        else:
            data[i, 1] = -1
            i+=1
    return data

def computeHist(filepaths, tonic=None, bins=None, normed=True, folded=False, vocal=True, weight="duration", quantized=False, slopeThresh=1500, centsThresh=50, durationThresh=None, annotationFile="/media/CompMusic/audio/users/gkoduri/Workspace/annotations.yaml"): 
    """Computes histogram of given array of filepaths. If the filepaths
    has more than one path, the returned histogram is an averaged one."""
    allCents = []
    for filepath in filepaths:
        print filepath
        mbid = filepath.split('/')[-1][:-4]
        if quantized:
            data = discreteFilter(mbid, srcdata=data, slopeThresh=slopeThresh, centsThresh=centsThresh, returnScale="cents", annotationFile=annotationFile, tonic=tonic)
            if durationThresh:
                data = durationFilter(data, durationThresh=durationThresh)
        else:
            #data = pitch(mbid, returnScale="cents", vocal=vocal, tonic=tonic, annotationFile=annotationFile)
            data = np.loadtxt(filepath)

        if data == None:
            raise Exception("Pitch data not received.")
        if vocal:
            data = vocalFilter(data, mbid)
        
        #Histogram
        validPitch = data[:,1]
        validPitch = [i for i in validPitch if i > -10000]
        if folded:
            validPitch = map(lambda x: int(x%1200), validPitch)

        tonic = None
        allCents.extend(validPitch)
    
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
        _sum = sum(n[:,1])
        for i in n[:,1]:
            m.append(1.0*i/_sum)
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
        leftBoundMax = p-maxDistribThresh
        rightBoundMax = p+maxDistribThresh
        leftBound = leftBoundMax
        rightBound = rightBoundMax
        nearestValleyIndex = findNearestIndex(valleys[0], p)
        if p > valleys[0][nearestValleyIndex]:
            leftBound = valleys[0][nearestValleyIndex]
            if len(valleys[0][nearestValleyIndex+1:]) == 0:
                rightBound = p+maxDistribThresh
            else:
                offset = nearestValleyIndex+1
                nearestValleyIndex = findNearestIndex(valleys[0][offset:], p)
                rightBound = valleys[0][offset+nearestValleyIndex]
        else:
            rightBound = valleys[0][nearestValleyIndex]
            if len(valleys[0][:nearestValleyIndex]) == 0:
                leftBound = p-maxDistribThresh
            else:
                nearestValleyIndex = findNearestIndex(valleys[0][:nearestValleyIndex], p)
                leftBound = valleys[0][nearestValleyIndex]
        
        #Think in terms of x-axis.
        #leftBound should be at least minDistribThresh less than p, and at max maxDistribThresh less than p
        #rightBound should be at least minDistribThresh greater than p, and at max maxDistribThresh greater than p
        if leftBound < leftBoundMax:
            leftBound = leftBoundMax
        elif leftBound > p-minDistribThresh:
            leftBound = p-minDistribThresh

        if rightBound > rightBoundMax:
            rightBound = rightBoundMax
        elif rightBound < p+minDistribThresh:
            rightBound = p+minDistribThresh
            
        #Bounds should be at equal distance on either side. If they are not, make them.
        newThresh = 0
        if p-leftBound < rightBound-p:
            newThresh = p-leftBound
        else:
            newThresh = rightBound-p
        leftBound = p-newThresh
        rightBound = p+newThresh
        
        #extract the distribution and estimate the parameters
        distribution = cents[cents>=leftBound]
        distribution = distribution[distribution<=rightBound]
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

def computeParameters(mbids, smoothingFactor=11, histogramsPath="/media/CompMusic/audio/users/gkoduri/Workspace/intonationLib/data/method-1/vocal-histograms/", pitchPath="/media/CompMusic/audio/users/gkoduri/Workspace/features/pitch/"):
    """This is a wrapper function to bulk calculate parameters for a
    large number of files (run for one raaga at a time)"""

    allParams = {}
    for mbid in mbids:
        print mbid
        path = histogramsPath+mbid+".pickle"
        if not exists(path):
            [n, binCenters, cents] = computeHist([pitchPath+mbid+".txt"])
        else:
            [n, binCenters, cents] = pickle.load(file(path))
        n = gaussian_filter(n, smoothingFactor)
        peakInfo = peaks.peaks(n, binCenters, method="hybrid", window=100, peakAmpThresh=0.00005, valleyThresh=0.00003)
        params = characterizePeaks(peakInfo["peaks"], peakInfo["valleys"], cents, maxDistribThresh=50, minDistribThresh=20)
        allParams[mbid] = params
    return allParams

# Functions for mean-window approach

def isolateSwaras(data, windowSize=150, hopSize=30):
    """
    data: mx2 array with time and pitch value in cents as columns.
    windowSize: the size of window over which means are calculated in
    milliseconds.
    hopSize : hop size in milliseconds.

    returns a dictionary swaraDistributions with swaras as keys and an array of
    corresponding pitches as values.
    """
    windowSize = windowSize/1000.0
    hopSize = hopSize/1000.0
    exposure = int(windowSize/hopSize)
    boundary = windowSize-hopSize
    finaleInd = findNearestIndex(data[:, 0], data[-1, 0]-boundary)

    justIntonation = [['Sa', 1.0], ['R1',16.0/15.0], ['R2/G1',9.0/8.0],\
    ['G2/R3',6.0/5.0], ['G3',5.0/4.0], ['M1',4.0/3.0], ['M2',64.0/45.0],\
    ['P',3.0/2.0], ['D1',8.0/5.0], ['D2/N1',5.0/3.0],\
    ['D3/N2',16.0/9.0], ['N3',15.0/8.0]]
    swaras = ["Sa_", "R1_", "R2/G1_", "G2/R3_", "G3_", "M1_", "M2_", "P_", "D1_", "D2/N1_", "D3/N2_", "N3_", "Sa", "R1", "R2/G1", "G2/R3", "G3", "M1", "M2", "P", "D1", "D2/N1", "D3/N2", "N3", "Sa^", "R1^", "R2/G1^", "G2/R3^", "G3^", "M1^", "M2^", "P^", "D1^", "D2/N1^", "D3/N2^", "N3^"]
    #swaras = ["Sa_", "R1_", "R2/G1_", "G2/R3_", "G3_", "M1_", "M2_", "P_", "D1_", "D2/N1_", "D3/N2_", "N3_", "Sa", "R1", "R2/G1", "G2/R3", "G3", "M1", "M2", "P", "D1", "D2/N1", "D3/N2", "N3", "Sa^", "R1^", "R2/G1^", "G2/R3^", "G3^", "M1^", "M2^", "P^", "D1^", "D2/N1^", "D3/N2^", "N3^", "Sa^^", "R1^^", "R2/G1^^", "G2/R3^^", "G3^^", "M1^^", "M2^^", "P^^", "D1^^", "D2/N1^^", "D3/N2^^", "N3^^"]
    JICents = []
    for interval in justIntonation:
        p = int(1200.0*np.log2(interval[1]))
        JICents.append(p)
        JICents.append(p-1200)
        JICents.append(p+1200)
        #JICents.append(p+2400)
    JICents.sort()

    startInd = 0
    #HARDCODED
    interval = 0.00290254832393
    windowStep = windowSize/interval
    hopStep = hopSize/interval
    endInd = windowStep
    swaraDistributions = {}
    _means = []
    while endInd < finaleInd:
        temp = data[startInd:endInd, 1][data[startInd:endInd, 1]>-10000]
        _means.append(np.mean(temp))
        startInd = startInd+hopStep
        endInd = startInd+windowStep

    for i in xrange(exposure, len(_means)-exposure+1):
        _median = np.median(_means[i-exposure:i])
        if _median < -5000: continue
        ind = findNearestIndex(_median, JICents)
        sliceEnd = (i-exposure)*hopStep+windowStep
        sliceBegin = sliceEnd-hopStep
        #print sliceBegin, sliceEnd, JICents[ind]
        #newPitch[sliceBegin:sliceEnd] = JICents[ind]
        if swaras[ind] in swaraDistributions.keys():
            swaraDistributions[swaras[ind]] = np.append(swaraDistributions[swaras[ind]], data[sliceBegin:sliceEnd,1])
        else:
            swaraDistributions[swaras[ind]] = data[sliceBegin:sliceEnd,1]
        
    for swara in swaraDistributions.keys():
        swaraDistributions[swara] = swaraDistributions[swara][swaraDistributions[swara] > -10000]
    return swaraDistributions

#Plotting functions

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
    labelFontSize = 20
    for interval in JICents:
        p.axvline(x=interval[1], ls='-.', c='g', lw='1')
        p.text(interval[1]+3, spacing+4*max(n)/5, interval[0], fontsize=labelFontSize)
        if interval[1]-2400 >= min(binCenters):
            p.axvline(x=interval[1]-2400, ls='--', c='r', lw='0.5')
            p.text(interval[1]-2400+3, spacing+4*max(n)/5, interval[0], fontsize=labelFontSize)
        if interval[1]-1200 >= min(binCenters):
            p.axvline(x=interval[1]-1200, ls=':', c='b', lw='0.5')
            p.text(interval[1]-1200+3, spacing+4*max(n)/5, interval[0], fontsize=labelFontSize)
        if interval[1]+1200 <= max(binCenters):
            p.axvline(x=interval[1]+1200, ls=':', c='b', lw='0.5')
            p.text(interval[1]+1200+3, spacing+4*max(n)/5, interval[0], fontsize=labelFontSize)
        if interval[1]+2400 <= max(binCenters):
            p.axvline(x=interval[1]+2400, ls='-.', c='r', lw='0.5')
            p.text(interval[1]+2400+3, spacing+4*max(n)/5, interval[0], fontsize=labelFontSize)
        spacing = -1*spacing;

def plotHist(n, binCenters, peakInfo=None, newFigure=True, intervals=True):
    """This function plots histogram together with its smoothed
    version and peak information if provided. Just intonation 
    intervals are plotted for a reference."""
    if newFigure:
        p.figure()
    #p.plot(binCenters, n, ls='-', c='#00CC00', lw='1')
    p.plot(binCenters, n, ls='-', c='b', lw='1.5')
    
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
    if intervals:
        labelFontSize = 15
        spacing = 0.02*max(n)
        for interval in JICents:
            p.axvline(x=interval[1], ls='-.', c='g', lw='1.5')
            #p.text(interval[1]+3, spacing+4*max(n)/5, interval[0], fontsize=labelFontSize)
            if interval[1]-2400 >= min(binCenters):
                p.axvline(x=interval[1]-2400, ls='--', c='r', lw='0.5')
                #p.text(interval[1]-2400+3, spacing+4*max(n)/5, interval[0], fontsize=labelFontSize)
            if interval[1]-1200 >= min(binCenters):
                p.axvline(x=interval[1]-1200, ls=':', c='b', lw='0.5')
                #p.text(interval[1]-1200+3, spacing+4*max(n)/5, interval[0], fontsize=labelFontSize)
            if interval[1]+1200 <= max(binCenters):
                p.axvline(x=interval[1]+1200, ls=':', c='b', lw='0.5')
                #p.text(interval[1]+1200+3, spacing+4*max(n)/5, interval[0], fontsize=labelFontSize)
            if interval[1]+2400 <= max(binCenters):
                p.axvline(x=interval[1]+2400, ls='-.', c='r', lw='0.5')
                #p.text(interval[1]+2400+3, spacing+4*max(n)/5, interval[0], fontsize=labelFontSize)
            spacing = -1*spacing;

    p.title("Tonic-aligned pitch histogram without folding octaves")
    p.xlabel("Pitch value (Cents)")
    p.ylabel("Normalized frequency of occurence")
    if peakInfo:
        p.plot(peakInfo["peaks"][0], peakInfo["peaks"][1], 'rD', ms=10)
        #p.plot(peakInfo["valleys"][0], peakInfo["valleys"][1], 'yD', ms=5)
    p.show()


#Functions to build database

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
        filewriter.write(str(raagas[i])+"\n")
    filewriter.close()
    return raagas

def filepathByMBID(mbid, pathfile="/media/CompMusic/Carnatic/metadata/Carnatic.yaml"):
    """
    filepathByMBID(mbid, pathfile="/media/CompMusic/Carnatic/metadata/Carnatic.yaml")

    Returns filepath of recording given by the mbid.
    """
    data = yaml.load(file(pathfile))
    if mbid in data.keys():
        return data[mbid]["path"]
    else:
        return None

def recordingsByRaagaAliases(raaga, sim_thresh=0.7, recursion=2):
    db = dao.DAO()
    tags = db.getTagsByCategory("raaga")
    raagas = [i["tag"] for i in tags[1]]
    aliases = stringDuplicates(raaga, raagas, simThresh=sim_thresh, recursion=recursion)
    recordings = {}
    #recordings[mbid] = path

    for alias in aliases:
        data = db.getRecordingsByTag("tag", alias, category = "raaga")
        if data[0]:
            for i in xrange(len(data[1])):
                uuid = data[1][i]['uuid']
                #path = filepathByMBID(uuid)
                recordings[uuid] = [alias]
    return recordings

