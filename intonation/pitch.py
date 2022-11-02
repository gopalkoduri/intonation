from __future__ import division
import numpy as np
from scipy.ndimage import median_filter
from . import utils


class Pitch:
    def __init__(self, timestamps, pitch):
        self.timestamps = timestamps
        self.pitch_raw = pitch
        self.pitch = pitch

    def reset(self):
        self.pitch = self.pitch_raw

    def discretize(self, intervals, slope_thresh=1500, cents_thresh=50):
        """
        This function takes the pitch data and returns it quantized to given
        set of intervals. All transactions must happen in cent scale.

        slope_thresh is the bound beyond which the pitch contour is said to transit
        from one svara to another. It is specified in cents/sec.

        cents_thresh is a limit within which two pitch values are considered the same.
        This is what pushes the quantization limit.

        The function returns quantized pitch data.
        """

        #eps = np.finfo(float).eps
        #pitch = median_filter(pitch, 7)+eps

        self.pitch = median_filter(self.pitch, 7)
        pitch_quantized = np.zeros(len(self.pitch))
        pitch_quantized[0] = utils.find_nearest_index(intervals, self.pitch[0])
        pitch_quantized[-1] = utils.find_nearest_index(intervals, self.pitch[-1])

        for i in xrange(1, len(self.pitch)-1):
            if self.pitch[i] == -10000:
                pitch_quantized[i] = -10000
                continue
            slope_back = abs((self.pitch[i] - self.pitch[i-1])/(self.timestamps[i] - self.timestamps[i-1]))
            slope_front = abs((self.pitch[i+1] - self.pitch[i])/(self.timestamps[i+1] - self.timestamps[i]))
            if slope_front < slope_thresh or slope_back < slope_thresh:
                ind = utils.find_nearest_index(intervals, self.pitch[i])
                cents_diff = abs(self.pitch[i] - intervals[ind])
                if cents_diff <= cents_thresh:
                    pitch_quantized[i] = intervals[ind]
                else:
                    pitch_quantized[i] = -10000
            else:
                pitch_quantized[i] = -10000

        self.pitch = pitch_quantized

    def enforce_duration(self, duration_thresh):
        """
        This method takes a quantized pitch contour and filters out
        those time sections where the contour is not long enough, as specified
        by duration threshold (given in milliseconds).

        All transactions assume data in cent scale.
        """
        i = 1
        while i < len(self.pitch)-1:
            if self.pitch[i] == -10000:
                i += 1
                continue
            if self.pitch[i]-self.pitch[i-1] != 0 and self.pitch[i+1]-self.pitch[i] == 0:
                start = i
                while i < len(self.pitch) and self.pitch[i+1]-self.pitch[i] == 0:
                    i += 1
                if (self.timestamps[i]-self.timestamps[start])*1000 < duration_thresh:
                    self.pitch[start:i+1] = np.zeros(i+1-start)-10000
            else:
                self.pitch[i] = -10000
                i += 1

    def fit_lines(self, window=1500, break_thresh=1500):
        """
        Fits lines to pitch contours.

        :param window: size of each chunk to which linear equation is to be fit (in milliseconds).
        To keep it simple, hop is chosen to be one third of the window.
        :param break_thresh: If there is silence beyond this limit (in milliseconds),
        the contour will be broken there into two so that we don't fit a line over and
        including the silent region.
        """
        window /= 1000
        hop = window/3
        break_thresh /= 1000

        #cut the whole song into pieces if there are gaps more than break_thresh seconds
        i = 0
        break_indices = []
        count = 0
        while i < len(self.pitch):
            if self.pitch[i] == -10000:
                count = 1
                start_index = i
                while i < len(self.pitch) and self.pitch[i] == -10000:
                    count += 1
                    i += 1
                end_index = i-1
                if self.timestamps[end_index]-self.timestamps[start_index] >= break_thresh:
                    break_indices.append([start_index, end_index])
            i += 1
        break_indices = np.array(break_indices)

        #In creating the data blocks which are not silences, note that we
        # take complimentary break indices. i.e., if [[s1, e1], [s2, e2] ...]
        # is break_indices, we take e1-s2, e2-s3 chunks and build data blocks

        data_blocks = []
        if len(break_indices) == 0:
            t_pitch = self.pitch.reshape(len(self.pitch), 1)
            t_timestamps = self.timestamps.reshape(len(self.timestamps), 1)
            data_blocks = [np.append(t_timestamps, t_pitch, axis=1)]
        else:
            if break_indices[0, 0] != 0:
                t_pitch = self.pitch[:break_indices[0, 0]]
                t_pitch = t_pitch.reshape(len(t_pitch), 1)
                t_timestamps = self.timestamps[:break_indices[0, 0]]
                t_timestamps = t_timestamps.reshape(len(t_timestamps), 1)
                data_blocks.append(np.append(t_timestamps, t_pitch, axis=1))
            block_start = break_indices[0, 1]
            for i in xrange(1, len(break_indices)):
                block_end = break_indices[i, 0]
                t_pitch = self.pitch[block_start:block_end]
                t_pitch = t_pitch.reshape(len(t_pitch), 1)
                t_timestamps = self.timestamps[block_start:block_end]
                t_timestamps = t_timestamps.reshape(len(t_timestamps), 1)
                data_blocks.append(np.append(t_timestamps, t_pitch, axis=1))
                block_start = break_indices[i, 1]
            if block_start != len(self.pitch)-1:
                t_pitch = self.pitch[block_start:]
                t_pitch = t_pitch.reshape(len(t_pitch), 1)
                t_timestamps = self.timestamps[block_start:]
                t_timestamps = t_timestamps.reshape(len(t_timestamps), 1)
                data_blocks.append(np.append(t_timestamps, t_pitch, axis=1))

        label_start_offset = (window-hop)/2
        label_end_offset = label_start_offset+hop

        #dataNew = np.zeros_like(data)
        #dataNew[:, 0] = data[:, 0]
        data_new = np.array([[0, 0]])
        for data in data_blocks:
            start_index = 0
            while start_index < len(data)-1:
                end_index = utils.find_nearest_index(data[:, 0], data[start_index][0]+window)
                segment = data[start_index:end_index]
                if len(segment) == 0:
                    start_index = utils.find_nearest_index(data[:, 0], data[start_index, 0]+hop)
                    continue
                segment_clean = np.delete(segment, np.where(segment[:, 1] == -10000), axis=0)
                if len(segment_clean) == 0:
                    #After splitting into blocks, this loop better not come into play
                    #raise ValueError("This part of the block is absolute silence! Make sure block_thresh >= window!")
                    start_index = utils.find_nearest_index(data[:, 0], data[start_index, 0]+hop)
                    continue
                n_clean = len(segment_clean)
                x_clean = np.matrix(segment_clean[:, 0]).reshape(n_clean, 1)
                y_clean = np.matrix(segment_clean[:, 1]).reshape(n_clean, 1)
                #return [x_clean, y_clean]
                theta = utils.normal_equation(x_clean, y_clean)

                #determine the start and end of the segment to be labelled
                label_start_index = utils.find_nearest_index(x_clean, data[start_index, 0]+label_start_offset)
                label_end_index = utils.find_nearest_index(x_clean, data[start_index, 0]+label_end_offset)
                x_clean = x_clean[label_start_index:label_end_index]
                #return x_clean
                x_clean = np.insert(x_clean, 0, np.ones(len(x_clean)), axis=1)
                newy = x_clean*theta
                result = np.append(x_clean[:, 1], newy, axis=1)
                data_new = np.append(data_new, result, axis=0)

                start_index = utils.find_nearest_index(data[:, 0], data[start_index, 0]+hop)

        return [data_new[:, 0], data_new[:, 1]]
