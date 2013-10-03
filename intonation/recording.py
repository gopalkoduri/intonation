from __future__ import division
import pickle
from pypeaks import Data
from scipy.stats import variation, skew, kurtosis
from pitch import Pitch
import utils
import numpy as np
#TODO: Vocal filter


class Recording:
    def __init__(self, pitch_obj):
        assert isinstance(self.pitch_obj, Pitch)
        self.pitch_obj = pitch_obj
        self.histogram = None
        self.intonation_profile = None
        self.contour_labels = None

    def serialize_hist(self, path):
        pickle.dump(self.histogram, file(path, 'w'))

    def serialize_intonation(self, path):
        pickle.dump(self.intonation_profile, file(path, 'w'))

    def serialize_contour_labels(self, path):
        pickle.dump(self.contour_labels, file(path, 'w'))

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
        :param weight: It can be one of the 'duration' or 'instance'. In the latter case, make
        sure that the pitch object has the pitch values discretized, and the
        intervals argument must be passed.
        :param intervals: an array of intervals.
        """
        #Step 1: get the right pitch values
        assert isinstance(self.pitch_obj.pitch, np.ndarray)
        valid_pitch = self.pitch_obj.pitch
        valid_pitch = [i for i in valid_pitch if i > -10000]
        if folded:
            valid_pitch = map(lambda x: int(x % 1200), valid_pitch)

        #Step 2: based on the weighing scheme, compute the histogram
        if weight == "duration":
            #Step 2.1 set the number of bins (if not passed)
            if not bins:
                bins = max(valid_pitch) - min(valid_pitch)
            n, bin_edges = np.histogram(valid_pitch, bins, density=density)
            bin_centers = 0.5 * (bin_edges[1:] + bin_edges[:-1])
            self.histogram = Data(bin_centers, n)
        elif weight == "instance":
            n = {}
            i = 1
            while i < len(valid_pitch) - 1:
                if (valid_pitch[i] - valid_pitch[i - 1] != 0) and \
                        (valid_pitch[i + 1] - valid_pitch[i] == 0):
                    #TODO: Do we need the intervals?? Remove them arguments if not necessary
                    #TODO: the following code needs to be removed if discritization works
                    # properly. If it doesn't there will be values like 1198.5, 1200, 1203.8
                    # etc.. for each interval, in which case we need this block.
                    #flag = 0
                    #for val in xrange(valid_pitch[i] - 3, valid_pitch[i] + 3, 1):
                    #    if val in n.keys():
                    #        n[val] += 1
                    #        flag = 1
                    #        break
                    #if flag == 0:
                    n[valid_pitch[i]] = 1
                i += 1
            n = n.items()
            n.sort(key=lambda x: x[0])
            n = np.array(n)

            median_diff = np.median(np.diff(n[:, 1]))
            bin_edges = [n[0, 1] - median_diff/2]
            bin_edges.extend(median_diff/2 + n[:, 1])
            n[:, 1] = n[:, 1]/(n[:, 1].sum()*np.diff(bin_edges))
            self.histogram = Data(n[:, 0], n[:, 1])

    def parametrize_peaks(self, intervals, max_peakwidth=50, min_peakwidth=25, symmetric_bounds=True):
        """
        Computes and stores the intonation profile of an audio recording.

        :param intervals: these will be the reference set of intervals to which peak positions
         correspond to. For each interval, the properties of corresponding peak, if exists,
         will be computed and stored as intonation profile.
        :param max_peakwidth: the maximum allowed width of the peak at the base for computing
        parameters of the distribution.
        :param min_peakwidth: the minimum allowed width of the peak at the base for computing
        parameters of the distribution.
        """
        assert isinstance(self.pitch_obj.pitch, np.ndarray)
        valid_pitch = self.pitch_obj.pitch
        valid_pitch = [i for i in valid_pitch if i > -10000]

        parameters = {}
        for i in xrange(len(self.histogram.peaks["peaks"][0])):
            peak_pos = self.histogram.peaks["peaks"][0][i]
            #Set left and right bounds of the distribution.
            max_leftbound = peak_pos - max_peakwidth
            max_rightbound = peak_pos + max_peakwidth
            leftbound = max_leftbound
            rightbound = max_rightbound
            nearest_valleyindex = utils.find_nearest_index(self.histogram.peaks["valleys"][0], peak_pos)
            if peak_pos > self.histogram.peaks["valleys"][0][nearest_valleyindex]:
                leftbound = self.histogram.peaks["valleys"][0][nearest_valleyindex]
                if len(self.histogram.peaks["valleys"][0][nearest_valleyindex + 1:]) == 0:
                    rightbound = peak_pos + max_peakwidth
                else:
                    offset = nearest_valleyindex + 1
                    nearest_valleyindex = utils.find_nearest_index(
                        self.histogram.peaks["valleys"][0][offset:], peak_pos)
                    rightbound = self.histogram.peaks["valleys"][0][offset + nearest_valleyindex]
            else:
                rightbound = self.histogram.peaks["valleys"][0][nearest_valleyindex]
                if len(self.histogram.peaks["valleys"][0][:nearest_valleyindex]) == 0:
                    leftbound = peak_pos - max_peakwidth
                else:
                    nearest_valleyindex = utils.find_nearest_index(
                        self.histogram.peaks["valleys"][0][:nearest_valleyindex], peak_pos)
                    leftbound = self.histogram.peaks["valleys"][0][nearest_valleyindex]

            #In terms of x-axis, leftbound should be at least min_peakwidth
            # less than peak_pos, and at max max_peakwidth less than peak_pos,
            # and viceversa for the rightbound.
            if leftbound < max_leftbound:
                leftbound = max_leftbound
            elif leftbound > peak_pos - min_peakwidth:
                leftbound = peak_pos - min_peakwidth

            if rightbound > max_rightbound:
                rightbound = max_rightbound
            elif rightbound < peak_pos + min_peakwidth:
                rightbound = peak_pos + min_peakwidth

            #If symmetric bounds are asked for, then make the bounds symmetric
            if symmetric_bounds:
                if peak_pos - leftbound < rightbound - peak_pos:
                    imbalance = (rightbound - peak_pos) - (peak_pos - leftbound)
                    rightbound -= imbalance
                else:
                    imbalance = (peak_pos - leftbound) - (rightbound - peak_pos)
                    leftbound += imbalance

            #extract the distribution and estimate the parameters
            distribution = valid_pitch[valid_pitch >= leftbound]
            distribution = distribution[distribution <= rightbound]
            print peak_pos, "\t", len(distribution), "\t", leftbound, "\t", rightbound

            interval_index = utils.find_nearest_index(intervals, peak_pos)
            interval = intervals[interval_index]
            _mean = float(np.mean(distribution))
            _variance = float(variation(distribution))
            _skew = float(skew(distribution))
            _kurtosis = float(kurtosis(distribution))
            pearson_skew = float(3.0 * (_mean - peak_pos) / np.sqrt(abs(_variance)))
            parameters[interval] = {"position": float(peak_pos),
                                    "mean": _mean,
                                    "amplitude": float(self.histogram.peaks[1][i]),
                                    "variance": _variance,
                                    "skew1": _skew,
                                    "skew2": pearson_skew,
                                    "kurtosis": _kurtosis}

        self.intonation_profile = parameters

    def label_contours(self, intervals, window=150, hop=30):
        """
        In a very flowy contour, it is not trivial to say which pitch value corresponds
         to what interval. This function labels pitch contours with intervals by guessing
         from the characteristics of the contour and its melodic context.

        :param window: the size of window over which the context is gauged, in milliseconds.
        :param hop: hop size in milliseconds.
        """
        window /= 1000.0
        hop /= 1000.0
        exposure = int(window / hop)

        boundary = window - hop
        final_index = utils.find_nearest_index(self.pitch_obj.timestamps,
                                               self.pitch_obj.timestamps[-1] - boundary)

        interval = np.median(np.diff(self.pitch_obj.timestamps))
        #interval = 0.00290254832393
        window_step = window / interval
        hop_step = hop / interval
        start_index = 0
        end_index = window_step
        contour_labels = {}
        means = []
        while end_index < final_index:
            temp = self.pitch_obj.pitch[start_index:end_index][self.pitch_obj.pitch[start_index:end_index] > -10000]
            means.append(np.mean(temp))
            start_index = start_index + hop_step
            end_index = start_index + window_step

        for i in xrange(exposure, len(means) - exposure + 1):
            _median = np.median(means[i - exposure:i])
            if _median < -5000:
                continue
            ind = utils.find_nearest_index(_median, intervals)
            contour_start = (i - exposure) * hop_step + window_step
            contour_end = contour_start - hop_step
            #print sliceBegin, sliceEnd, JICents[ind]
            #newPitch[sliceBegin:sliceEnd] = JICents[ind]
            if intervals[ind] in contour_labels.keys():
                contour_labels[intervals[ind]].append([contour_start, contour_end])
            else:
                contour_labels[intervals[ind]] = [[contour_start, contour_end]]

        self.contour_labels = contour_labels