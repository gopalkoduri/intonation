import numpy as np
from scipy.ndimage import median_filter
import utils


class Pitch:
    def __init__(self, timestamps, pitch):
        self.timestamps = timestamps
        self.pitch = pitch

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
