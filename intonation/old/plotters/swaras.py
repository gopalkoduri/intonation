#!/usr/bin/env python

import pickle
import yaml
import sys
import numpy as np
from scipy.ndimage.filters import gaussian_filter
from pylab import *

raagas = ["abhogi", "begada", "mohanam"]
raagaMBIDs = yaml.load(file(sys.argv[1]))

f, ax = subplots(3, sharex=True)

styles = ["-", "--", "-", "--", "-", "--"]
colors = ["#000000", "#550000", "#009900", "#0000BB", "#CCCC00", "#996600"]

raagaCount = 0
recCount = 0
for r in raagas:
    recCount = 0
    for mbid in raagaMBIDs[r][:5]:
        data = pickle.load(file(sys.argv[2]+"/"+mbid+".pickle"))
        [n, be] = np.histogram(data["G3"], bins=1200, density=True)
        bc = (be[:-1]+be[1:])/2.0
        ns = gaussian_filter(n, 7)
        ax[raagaCount].plot(bc, ns, color=colors[recCount], ls=styles[recCount], lw="1.25")
        recCount += 1
        
    ax[raagaCount].grid(color="0.35")
    raagaCount += 1
