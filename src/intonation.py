from pypeaks import Data, Intervals
from pitch import Pitch

import numpy as np

#TODO: Vocal filter

class Intonation:
    def __init__(self, pitch_obj):
        assert isinstance(self.pitch_obj, Pitch)
        self.pitch_obj = pitch_obj
        self.histogram = None

    def compute_hist(self, bins=None, density=True, folded=False, weight="duration",
                     intervals=None):
        """
        Computes histogram from the pitch data in Pitch object (pitch), and creates
        a Data object (pypeaks).

        :param bins: Refers to number of bins in the histogram, determines the granularity.
        If it is not set, the number of bins which gives the highest granularity is chosen
        automatically.
        :param density: defaults to True, which means the histogram will be a normalized one.
        :param folded: defaults to False. When set to True, all the octaves are folded to one.
        :param weight: It can be one of the 'duration' or 'instance'. In the latter case, the
        number of bins must be stated, and they should match with the Interval object passed.
        """
        #Step 1: get the right pitch values
        assert isinstance(self.pitch_obj.pitch, np.ndarray)
        valid_pitch = self.pitch_obj.pitch
        valid_pitch = [i for i in valid_pitch if i > -10000]
        if folded:
            valid_pitch = map(lambda x: int(x % 1200), valid_pitch)

        #Step 2: set the number of bins (if not passed)
        if not bins:
            bins = max(valid_pitch) - min(valid_pitch)

        #Step 3: based on the weighing scheme, compute the histogram
        if weight == "duration":
            n, bin_edges = np.histogram(valid_pitch, bins, density=density)
            bin_centers = 0.5 * (bin_edges[1:] + bin_edges[:-1])
            
        elif weight == "instance":
            if bins != 12:
                print "Quantization with bins != 12 is not handled. Quitting."
                exit()
            n = {}
            i = 1
            while (i < len(valid_pitch) - 1):
                if (valid_pitch[i] - valid_pitch[i - 1] != 0) and (valid_pitch[i + 1] - valid_pitch[i] == 0):
                    flag = 0
                    for val in xrange(valid_pitch[i] - 3, valid_pitch[i] + 3, 1):
                        if val in n.keys():
                            n[val] += 1
                            flag = 1
                            break
                    if flag == 0:
                        n[valid_pitch[i]] = 1
                i += 1
            n = n.items()
            n.sort(key=lambda x: x[0])
            n = np.array(n)
            m = []
            _sum = sum(n[:, 1])
            for i in n[:, 1]:
                m.append(1.0 * i / _sum)
            return [np.array(m), n[:, 0], valid_pitch]

    def characterizePeaks(peaks, valleys, cents, maxDistribThresh=50, minDistribThresh=25):
        """Given peaks and cent values, it returns the characteristics
        of each peak."""
        justIntonation = [['Sa', 1.0], ['R1', 16.0 / 15.0], ['R2/G1', 9.0 / 8.0], \
                          ['G2/R3', 6.0 / 5.0], ['G3', 5.0 / 4.0], ['M1', 4.0 / 3.0], ['M2', 64.0 / 45.0], \
                          ['P', 3.0 / 2.0], ['D1', 8.0 / 5.0], ['D2/N1', 5.0 / 3.0], \
                          ['D3/N2', 16.0 / 9.0], ['N3', 15.0 / 8.0]]
        swaras = ["Sa_", "R1_", "R2/G1_", "G2/R3_", "G3_", "M1_", "M2_", "P_", "D1_", "D2/N1_", "D3/N2_", "N3_", "Sa",
                  "R1", "R2/G1", "G2/R3", "G3", "M1", "M2", "P", "D1", "D2/N1", "D3/N2", "N3", "Sa^", "R1^", "R2/G1^",
                  "G2/R3^", "G3^", "M1^", "M2^", "P^", "D1^", "D2/N1^", "D3/N2^", "N3^", "Sa^^", "R1^^", "R2/G1^^",
                  "G2/R3^^", "G3^^", "M1^^", "M2^^", "P^^", "D1^^", "D2/N1^^", "D3/N2^^", "N3^^"]
        JICents = []
        for interval in justIntonation:
            p = round(1200.0 * np.log2(interval[1]), 2)
            JICents.append(p)
            JICents.append(p - 1200)
            JICents.append(p + 1200)
            JICents.append(p + 2400)

        JICents.sort()
        cents = np.array(cents)
        parameters = {}
        for i in xrange(len(peaks[0])):
            p = peaks[0][i]
            #Set left and right bounds of the distribution.
            leftBoundMax = p - maxDistribThresh
            rightBoundMax = p + maxDistribThresh
            leftBound = leftBoundMax
            rightBound = rightBoundMax
            nearestValleyIndex = find_nearest_index(valleys[0], p)
            if p > valleys[0][nearestValleyIndex]:
                leftBound = valleys[0][nearestValleyIndex]
                if len(valleys[0][nearestValleyIndex + 1:]) == 0:
                    rightBound = p + maxDistribThresh
                else:
                    offset = nearestValleyIndex + 1
                    nearestValleyIndex = find_nearest_index(valleys[0][offset:], p)
                    rightBound = valleys[0][offset + nearestValleyIndex]
            else:
                rightBound = valleys[0][nearestValleyIndex]
                if len(valleys[0][:nearestValleyIndex]) == 0:
                    leftBound = p - maxDistribThresh
                else:
                    nearestValleyIndex = find_nearest_index(valleys[0][:nearestValleyIndex], p)
                    leftBound = valleys[0][nearestValleyIndex]

            #Think in terms of x-axis.
            #leftBound should be at least minDistribThresh less than p, and at max maxDistribThresh less than p
            #rightBound should be at least minDistribThresh greater than p, and at max maxDistribThresh greater than p
            if leftBound < leftBoundMax:
                leftBound = leftBoundMax
            elif leftBound > p - minDistribThresh:
                leftBound = p - minDistribThresh

            if rightBound > rightBoundMax:
                rightBound = rightBoundMax
            elif rightBound < p + minDistribThresh:
                rightBound = p + minDistribThresh

            #Bounds should be at equal distance on either side. If they are not, make them.
            newThresh = 0
            if p - leftBound < rightBound - p:
                newThresh = p - leftBound
            else:
                newThresh = rightBound - p
            leftBound = p - newThresh
            rightBound = p + newThresh

            #extract the distribution and estimate the parameters
            distribution = cents[cents >= leftBound]
            distribution = distribution[distribution <= rightBound]
            #print p, "\t", len(distribution), "\t", leftBound, "\t", rightBound
            swaraIndex = find_nearest_index(JICents, p)
            swara = swaras[swaraIndex]
            _mean = float(np.mean(distribution))
            _variance = float(variation(distribution))
            _skew = float(skew(distribution))
            _kurtosis = float(kurtosis(distribution))
            pearsonSkew = float(3.0 * (_mean - p) / np.sqrt(abs(_variance)))
            parameters[swara] = {"position": float(p), "mean": _mean, "variance": _variance, "skew1": _skew,
                                 "kurtosis": _kurtosis, "amplitude": float(peaks[1][i]), "skew2": pearsonSkew}
        return parameters

    def computeParameters(mbids, smoothingFactor=11,
                          histogramsPath="/media/CompMusic/audio/users/gkoduri/Workspace/intonationLib/data/method-1/vocal-histograms/",
                          pitchPath="/media/CompMusic/audio/users/gkoduri/Workspace/features/pitch/"):
        """This is a wrapper function to bulk calculate parameters for a
        large number of files (run for one raaga at a time)"""

        allParams = {}
        for mbid in mbids:
            print mbid
            path = histogramsPath + mbid + ".pickle"
            if not exists(path):
                [n, binCenters, cents] = computeHist([pitchPath + mbid + ".txt"])
            else:
                [n, binCenters, cents] = pickle.load(file(path))
            n = gaussian_filter(n, smoothingFactor)
            peakInfo = peaks.peaks(n, binCenters, method="hybrid", window=100, peakAmpThresh=0.00005,
                                   valleyThresh=0.00003)
            params = characterizePeaks(peakInfo["peaks"], peakInfo["valleys"], cents, maxDistribThresh=50,
                                       minDistribThresh=20)
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
        windowSize = windowSize / 1000.0
        hopSize = hopSize / 1000.0
        exposure = int(windowSize / hopSize)
        boundary = windowSize - hopSize
        finaleInd = find_nearest_index(data[:, 0], data[-1, 0] - boundary)

        justIntonation = [['Sa', 1.0], ['R1', 16.0 / 15.0], ['R2/G1', 9.0 / 8.0], \
                          ['G2/R3', 6.0 / 5.0], ['G3', 5.0 / 4.0], ['M1', 4.0 / 3.0], ['M2', 64.0 / 45.0], \
                          ['P', 3.0 / 2.0], ['D1', 8.0 / 5.0], ['D2/N1', 5.0 / 3.0], \
                          ['D3/N2', 16.0 / 9.0], ['N3', 15.0 / 8.0]]
        swaras = ["Sa_", "R1_", "R2/G1_", "G2/R3_", "G3_", "M1_", "M2_", "P_", "D1_", "D2/N1_", "D3/N2_", "N3_", "Sa",
                  "R1", "R2/G1", "G2/R3", "G3", "M1", "M2", "P", "D1", "D2/N1", "D3/N2", "N3", "Sa^", "R1^", "R2/G1^",
                  "G2/R3^", "G3^", "M1^", "M2^", "P^", "D1^", "D2/N1^", "D3/N2^", "N3^"]
        #swaras = ["Sa_", "R1_", "R2/G1_", "G2/R3_", "G3_", "M1_", "M2_", "P_", "D1_", "D2/N1_", "D3/N2_", "N3_", "Sa", "R1", "R2/G1", "G2/R3", "G3", "M1", "M2", "P", "D1", "D2/N1", "D3/N2", "N3", "Sa^", "R1^", "R2/G1^", "G2/R3^", "G3^", "M1^", "M2^", "P^", "D1^", "D2/N1^", "D3/N2^", "N3^", "Sa^^", "R1^^", "R2/G1^^", "G2/R3^^", "G3^^", "M1^^", "M2^^", "P^^", "D1^^", "D2/N1^^", "D3/N2^^", "N3^^"]
        JICents = []
        for interval in justIntonation:
            p = int(1200.0 * np.log2(interval[1]))
            JICents.append(p)
            JICents.append(p - 1200)
            JICents.append(p + 1200)
            #JICents.append(p+2400)
        JICents.sort()

        startInd = 0
        #HARDCODED
        interval = 0.00290254832393
        windowStep = windowSize / interval
        hopStep = hopSize / interval
        endInd = windowStep
        swaraDistributions = {}
        _means = []
        while endInd < finaleInd:
            temp = data[startInd:endInd, 1][data[startInd:endInd, 1] > -10000]
            _means.append(np.mean(temp))
            startInd = startInd + hopStep
            endInd = startInd + windowStep

        for i in xrange(exposure, len(_means) - exposure + 1):
            _median = np.median(_means[i - exposure:i])
            if _median < -5000: continue
            ind = find_nearest_index(_median, JICents)
            sliceEnd = (i - exposure) * hopStep + windowStep
            sliceBegin = sliceEnd - hopStep
            #print sliceBegin, sliceEnd, JICents[ind]
            #newPitch[sliceBegin:sliceEnd] = JICents[ind]
            if swaras[ind] in swaraDistributions.keys():
                swaraDistributions[swaras[ind]] = np.append(swaraDistributions[swaras[ind]],
                                                            data[sliceBegin:sliceEnd, 1])
            else:
                swaraDistributions[swaras[ind]] = data[sliceBegin:sliceEnd, 1]

        for swara in swaraDistributions.keys():
            swaraDistributions[swara] = swaraDistributions[swara][swaraDistributions[swara] > -10000]
        return swaraDistributions

        #Plotting functions
