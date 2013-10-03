# -*- coding:utf-8 -*-
#!/usr/bin/env python

import pickle
from scipy.ndimage.filters import gaussian_filter
from os import listdir
from pylab import *

#matplotlib.rcParams['ps.useafm'] = True
#matplotlib.rcParams['pdf.use14corefonts'] = True
#matplotlib.rcParams['text.usetex'] = True
#matplotlib.rcParams['text.latex.unicode']=True

#artists = ["vignesh", "ramakrishnamurthy", "dharini", "sreevidya", "prasanna"]
#raagas = ["abhogi"]

#artists = ["vignesh", "sreevidya", "prasanna"]
#raagas = ["begada"]

#artists = ["vignesh", "ramakrishnamurthy", "badrinarayanan", "prasanna"]
#raagas = ["kalyani"]
smoothingFactor = 5

#raagas = listdir("/home/gopal/workspace/intonationLib/data/varnam-analysis/recorded/distributions/")
#raagas = ['abhogi', 'begada', 'kalyani', 'mohanam', 'sahana', 'saveri', 'shree']
raagas = ['abhogi', 'begada', 'mohanam']
styles = ["-", "--", "-", "--", "-", "--"]
colors = ["#000000", "#550000", "#009900", "#0000BB", "#CCCC00", "#996600"]

distributions = {}
for r in raagas:
    artists =  listdir("/home/gopal/workspace/intonationLib/data/varnam-analysis/recorded/distributions/"+r)
    artists = [i[:-7] for i in artists]
    for a in artists:
        distributions[r+"_"+a] = pickle.load(file("/home/gopal/workspace/intonationLib/data/varnam-analysis/recorded/distributions/"+r+"/"+a+".pickle"))

for r in raagas:
    artists =  listdir("/home/gopal/workspace/intonationLib/data/varnam-analysis/recorded/distributions/"+r)
    artists = [i[:-7] for i in artists]
    artists.remove("dharini")
    for s in ["G"]:#distributions[r+"_"+artists[0]].keys():
        figure(figsize=(16, 9))
        #xlabel(r.title()+" : "+s)
        legendMarkers = []
        count = 0
        for a in artists:
            legendMarkers.append(a)
            [n, binEdges] = histogram(distributions[r+"_"+a][s], bins=1200, density=True)
            nS = gaussian_filter(n, 5)
            binCenters = 0.5*(binEdges[1:]+binEdges[:-1])
            plot(binCenters, nS, color=colors[count], ls=styles[count], lw="1.25")
            #plot(binCenters, nS)
            count += 1

        _min = (int(min(binCenters))/100-1)*100
        _max = (int(max(binCenters))/100+2)*100
        _xticks = array(xrange(_min, _max, 100))
        _xticks = _xticks.astype("int")
        xticks(_xticks, fontsize=24)
        #legend(legendMarkers, fontsize=24)
        yticks(fontsize=24)
        xlim(-100, 1300)
        ylim(ymin=0.00001)
        grid(color="0.35")
        savefig("/home/gopal/workspace/intonationLib/data/varnam-analysis/recorded/plots/"+r+"-"+s+".pdf",\
               bbox_inches=0)
