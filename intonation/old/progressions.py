from __future__ import division
from scipy.stats import linregress
import numpy as np
import linearRegression as lR
import utils

def fitLines(data, window=1500, hop=500, breakThresh=1500):
    """
    mbid: MBID of the recording
    window: Size of each chunk to which linear eqn is to be fit (in millisec)
    """

    window = window/1000
    hop = hop/1000
    breakThresh = breakThresh/1000

    labelOffsetStart = (window-hop)/2
    labelOffsetEnd = labelOffsetStart+hop

    numSamples = len(data)
    #data = np.delete(data, np.where(data[:, 1] == -10000), axis=0)
    
    #cut the whole song into pieces given there are gaps more than breakThresh seconds
    i = 0
    breakIndices = []
    while i < numSamples:
        if data[i, 1] == -10000:
            count = 1
            startInd = i
            while (i < numSamples) and (data[i, 1] == -10000):
                count+=1
                i+=1
            endInd = i-1
            if data[endInd, 0]-data[startInd, 0] >= breakThresh:
                breakIndices.append([startInd, endInd])
        i+=1
    breakIndices = np.array(breakIndices)

    dataBlocks = []
    if breakIndices[0, 0] != 0:
        dataBlocks.append(data[:breakIndices[0, 0]])
    blockStart = breakIndices[0, 1]
    for i in xrange(1, len(breakIndices)):
        blockEnd = breakIndices[i, 0]
        dataBlocks.append(data[blockStart:blockEnd])
        blockStart = breakIndices[i, 1]
    if blockStart != numSamples-1:
        dataBlocks.append(data[blockStart:])

    #dataNew = np.zeros_like(data)
    #dataNew[:, 0] = data[:, 0]
    dataNew = np.array([[0, 0]])
    for data in dataBlocks:
        numSamples = len(data)
        startInd = 0
        while startInd < numSamples-1:
            endInd = utils.find_nearest_index(data[:,0], data[startInd][0]+window)

            segment = data[startInd:endInd]
            n = len(segment)
            segmentClean = np.delete(segment, np.where(segment[:, 1] == -10000), axis=0)
            if len(segmentClean) == 0:
                #After splitting into blocks, this loop should never come into
                #action
                print "Blocks does not seem to be well split!"
                startInd = utils.find_nearest_index(data[:, 0], data[startInd, 0]+hop)
                continue
            nClean = len(segmentClean)
            XClean = np.matrix(segmentClean[:, 0]).reshape(nClean, 1)
            yClean = np.matrix(segmentClean[:, 1]).reshape(nClean, 1)
            theta = lR.normalEquation(XClean, yClean)

            #determine the start and end of the segment to be labelled
            labelStartInd = utils.find_nearest_index(XClean,\
                                                data[startInd, 0]+labelOffsetStart)
            labelEndInd = utils.find_nearest_index(XClean,\
                                                data[startInd, 0]+labelOffsetEnd)
            XClean = XClean[labelStartInd:labelEndInd]
            nClean = len(XClean)
            X = np.insert(XClean, 0, np.ones(nClean), axis=1)

            XClean = XClean.reshape(1, XClean.size).A[0]
            #newy = X*theta
            #newy = newy.reshape(1, newy.size).A[0]
            newy = [theta[1]]*len(XClean)
            result = np.array([XClean, newy]).T
            dataNew = np.append(dataNew, result, axis=0)

            startInd = utils.find_nearest_index(data[:, 0], data[startInd, 0]+hop)

    return dataNew

