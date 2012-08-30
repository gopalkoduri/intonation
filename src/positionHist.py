from numpy import *
import yaml
from intonationLib import findNearestIndex

allParameters = yaml.load(file("../data/allParameters-25-50.yaml"))
param = "position"
positions = {}
allPositions = []
for raaga in allParameters.keys():
	positions[raaga] = []
	for rec in allParameters[raaga].keys():
		offset = 0
		#if "Sa^" in allParameters[raaga][rec].keys():
			#offset = allParameters[raaga][rec]["Sa^"][param]-1200
		for swara in allParameters[raaga][rec].keys():
			positions[raaga].append(allParameters[raaga][rec][swara][param]-offset)
			allPositions.append(allParameters[raaga][rec][swara][param]-offset)

#raaga = "khamas"
template = array(allPositions)
#template = array(positions[raaga])
template = template[template>-50]
template = template[template<2350]

equi = []
means = []
justRef = [0.0, 111.73, 203.91, 315.63999999999999, 386.31, 498.04000000000002, 609.77999999999997, 701.96000000000004, 813.69000000000005, 884.36000000000001, 996.09000000000003, 1088.27, 1200.0, 1311.73, 1403.9100000000001, 1515.6399999999999, 1586.3099999999999, 1698.04, 1809.78, 1901.96, 2013.6900000000001, 2084.3600000000001, 2196.0900000000001, 2288.27]
just = []
swaraRef = ['Sa', 'R1', 'R2/G1', 'G2/R3', 'G3', 'M1', 'M2', 'P', 'D1', 'D2/N1', 'D3/N2', 'N3', 'Sa^\wedge', 'R1^\wedge', 'R2/G1^\wedge', 'G2/R3^\wedge', 'G3^\wedge', 'M1^\wedge', 'M2^\wedge', 'P^\wedge', 'D1^\wedge', 'D2/N1^\wedge', 'D3/N2^\wedge', 'N3^\wedge']

for i in xrange(24):
	position = i*100
	peak = template[template>position-50]
	peak = peak[peak<position+50]
	if len(peak) > 20:
		#print position, "\t", len(peak), "\t", mean(peak), "\t", mode(peak), "\t", variation(peak), "\t", skew(peak)
		_mean = round(mean(peak), 2)
		De = 0
		Dj = 0
		for val in peak:
			De += abs(position-val)
			Dj += abs(justRef[findNearestIndex(justRef, position)]-val)
		De = round(De*1.0/len(peak), 2)
		Dj = round(Dj*1.0/len(peak), 2)
		#print "$"+swaraRef[findNearestIndex(justRef, position)]+"$", "&", _mean, "&", position-_mean, "&", justRef[findNearestIndex(justRef, position)]-_mean, "&", len(peak), "\\\\"
		print "$"+swaraRef[findNearestIndex(justRef, position)]+"$", "&", _mean, "&", De, "&", Dj, "&", len(peak), "\\\\"
		#print "\hline"
		equi.append(position)
		means.append(mean(peak))
		just.append(justRef[findNearestIndex(justRef, position)])