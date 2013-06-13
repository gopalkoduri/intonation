#!/usr/bin/env python

import essentia.standard as es
import numpy as np
import peakDetectionAlgos as pda

def findNearestIndex(arr,value):
    """For a given value, the function finds the nearest value 
    in the array and returns its index."""
    arr = np.array(arr)
    index=(np.abs(arr-value)).argmin()
    return index

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

    index = findNearestIndex(JICents, ET)
    return JICents[index]

def peaksBySlope(n, binCenters, lookahead=20, delta=0.00003, averageHist=False, refPeaks=None, refThresh=25):
    #In the average histogram, we get the refPeaks and refValleys which are then used to clean the peaks we get from a single performance. For average histogram, lookahead should be ~ 20, delta = 0.000025 to pick up only the prominent peaks. For histogram of a performance, lookahead should be ~=15, delta = 0.00003 to pick up even the little peaks. 
    _max, _min = pda.peakdetect(n, binCenters, lookahead, delta)
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
        propThresh = 30 #NOTE: Hardcoded
        tempPeaks = [i+1200 for i in refPeaks[0]]
        tempPeaks.extend([i-1200 for i in refPeaks[0]])
        extendedPeaks = []
        extendedPeaks.extend(refPeaks[0])
        for i in tempPeaks:
            #if a peak exists around, don't add this new one.
            nearestInd = findNearestIndex(refPeaks[0], i)
            diff = abs(refPeaks[0][nearestInd] - i)
            diff = np.mod(diff, 1200)
            if diff > propThresh:
                extendedPeaks.append(i)
        #print extendedPeaks
        for peakLocationIndex in xrange(len(xPeaks)):
            extPeakLocationIndex = findNearestIndex(extendedPeaks, xPeaks[peakLocationIndex])
            diff = abs(xPeaks[peakLocationIndex] - extendedPeaks[extPeakLocationIndex])
            diff = np.mod(diff, 1200)
            #print xPeaks[peakLocationIndex], extendedPeaks[extPeakLocationIndex], diff
            if diff < refThresh:
                xCleanPeaks.append(xPeaks[peakLocationIndex])
                yCleanPeaks.append(yPeaks[peakLocationIndex])
        return {"peaks":[xCleanPeaks, yCleanPeaks], "valleys":[xValleys, yValleys]}

def peaks(n, binCenters, method="JI", window=100, peakAmpThresh=0.00005, valleyThresh=0.00003):
    """
    This function expects smoothed histogram (i.e., n).

    method can be JI/ET/slope/hybrid.
    JI and ET methods do not use generic peak detection algorithm. They use intervals and window
    to pick up the local maximum, and later filter out irrelevant peaks by using empirically found
    thresholds. Slope approach uses generic peak picking algorithm which first finds peaks by slope
    and then applies bounds. Hybrid approach first finds peaks using generic peak detection algo, 
    then filters the peaks heuristically as in JI/ET.
    
    window refers to the cent range used while picking up the maxima.

    The method returns:
    {"peaks":[[peak positions], [peak amplitudes]], "valleys": [[valley positions], [valley amplitudes]]}
    """
    data = zip(binCenters, n)
    binCenters = np.array(binCenters)
    firstCenter = (min(binCenters)+1.5*window)/window*window
    lastCenter = (max(binCenters)-window)/window*window
    if firstCenter < -1200: firstCenter = -1200
    if lastCenter > 3600: lastCenter = 3600


    if method == "slope" or method == "hybrid":
        peaks = {}
        peakInfo = peaksBySlope(n, binCenters, lookahead=20, delta=valleyThresh, averageHist=True)

        #find correspondences between peaks and valleys, and set valleys are left and right Indices
        #see the other method(s) for clarity!

        peakData = peakInfo["peaks"]
        valleyData = peakInfo["valleys"]

        #print len(peakData[0]), len(peakData[1])
        for i in xrange(len(peakData[0])):
            nearestIndex = findNearestIndex(valleyData[0], peakData[0][i])
            if valleyData[0][nearestIndex] < peakData[0][i]:
                leftIndex = findNearestIndex(binCenters, valleyData[0][nearestIndex])
                if (len(valleyData[0][nearestIndex+1:]) == 0):
                    rightIndex = findNearestIndex(binCenters, peakData[0][i]+window/2.0)
                else:
                    offset = nearestIndex+1
                    nearestIndex = offset+findNearestIndex(valleyData[0][offset:], peakData[0][i])
                    rightIndex = findNearestIndex(binCenters, valleyData[0][nearestIndex])
            else:
                rightIndex = findNearestIndex(binCenters, valleyData[0][nearestIndex])
                if (len(valleyData[0][:nearestIndex]) == 0):
                    leftIndex = findNearestIndex(binCenters, peakData[0][i]-window/2.0)
                else:
                    nearestIndex = findNearestIndex(valleyData[0][:nearestIndex], peakData[0][i])
                    leftIndex = findNearestIndex(binCenters, valleyData[0][nearestIndex])

            pos = findNearestIndex(binCenters, peakData[0][i])
            #print binCenters[pos], peakData[1][i], binCenters[leftIndex], binCenters[rightIndex]
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
        while interval < lastCenter:
            if method == "ET":
                leftIndex = findNearestIndex(binCenters, interval-window/2)
                rightIndex = findNearestIndex(binCenters, interval+window/2)
                interval += window
            elif method == "JI" or method == "hybrid":
                leftIndex = findNearestIndex(binCenters, (interval+prevInterval)/2.0)
                prevInterval = interval
                interval = nextJI(interval)
                rightIndex = findNearestIndex(binCenters, (interval+prevInterval)/2.0)
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
            nearIndex = findNearestIndex(p2, p)
            if abs(p-p2[nearIndex]) < window/2.0: p2.pop(nearIndex)
    
        for p in p1: allPeaks[p] = slopePeaks[p]
        for p in p2: allPeaks[p] = peaks[p]
        peaks = allPeaks

    #Filter the peaks and retain eligible peaks, also get their valley points.

    # ----> peakAmpThresh <---- : remove the peaks which are below that

    for pos in peaks.keys():
        #pos is an index in binCenters/n. DOES NOT refer to a cent value.
        if peaks[pos][0] < peakAmpThresh:
            #print "peakAmp: ", binCenters[pos]
            peaks.pop(pos)

    #Check if either left or right valley is deeper than ----> valleyThresh <----.
    valleys = {}
    for pos in peaks.keys():
        leftLobe = n[peaks[pos][1]:pos]
        rightLobe = n[pos:peaks[pos][2]]
        #Sanity check: Is it a genuine peak? Size of distributions on either side of the peak should be comparable.
        if len(leftLobe) == 0 or len(rightLobe) == 0:
            continue
        if 1.0*len(leftLobe)/len(rightLobe) < 0.15 or 1.0*len(leftLobe)/len(rightLobe) > 6.67:
            #print "size: ", binCenters[pos]
            #peaks.pop(pos)
            continue

        leftValleyPos = np.argmin(leftLobe)
        rightValleyPos = np.argmin(rightLobe)
        if (abs(leftLobe[leftValleyPos]-n[pos]) < valleyThresh and abs(rightLobe[rightValleyPos]-n[pos]) < valleyThresh):
            #print "valley: ", binCenters[pos]
            peaks.pop(pos)
        else:
            valleys[peaks[pos][1]+leftValleyPos] = leftLobe[leftValleyPos]
            valleys[pos+rightValleyPos] = rightLobe[rightValleyPos]
    
    if len(peaks) > 0:
        temp1 = np.array(peaks.values())
        temp1 = temp1[:, 0]

        return {'peaks':[binCenters[peaks.keys()], temp1], 'valleys':[binCenters[valleys.keys()], valleys.values()]}
    else:
        return {'peaks':[[], []], 'valleys':[[], []]}

def extendPeaks(srcPeaks, propThresh=30):
    """Each peak in srcPeaks is checked for its presence in other octaves. If it does not exist, it is created. propThresh is the cent range within which the peak in the other octave is expected to be present, i.e., only if there is a peak within this cent range in other octaves, then the peak is considered to be present in that octave.
    """
    #octave propagation of the reference peaks
    tempPeaks = [i+1200 for i in srcPeaks["peaks"][0]]
    tempPeaks.extend([i-1200 for i in srcPeaks["peaks"][0]])
    extendedPeaks = []
    extendedPeaks.extend(srcPeaks["peaks"][0])
    for i in tempPeaks:
        #if a peak exists around, don't add this new one.
        nearestInd = findNearestIndex(srcPeaks["peaks"][0], i)
        diff = abs(srcPeaks["peaks"][0][nearestInd] - i)
        diff = np.mod(diff, 1200)
        if diff > propThresh:
            extendedPeaks.append(i)
    return extendedPeaks

