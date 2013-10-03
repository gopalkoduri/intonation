from numpy import loadtxt
import pitch as p
reload(p)
import recording as r
reload(r)

from pypeaks import Data, Intervals
import numpy as np

data = loadtxt("../../../data/features/pitch/88d8196a-123a-4306-9856-4ce3faca14fc.txt")
data = data[:, :2]

pobj = p.Pitch(data[:, 0], data[:, 1])
pobj.fit_lines()
