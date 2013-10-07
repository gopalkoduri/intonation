import pickle
from numpy import loadtxt
import pitch as p
reload(p)
import recording as r
reload(r)

from pypeaks import Data, Intervals
import numpy as np

data = loadtxt("examples/88d8196a-123a-4306-9856-4ce3faca14fc.txt")
data = data[:, :2]

intervals = pickle.load(file("examples/ji-intervals.pickle"))

pobj = p.Pitch(data[:, 0], data[:, 1])
rec = r.Recording(pobj)
rec.compute_hist(weight='duration')
rec.histogram.get_peaks()
rec.parametrize_peaks(intervals)
