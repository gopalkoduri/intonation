import yaml
import intonationLib as iL
from os.path import exists

data = yaml.load(file("pantuvarali.yaml"))
mbids = data.keys()
filepaths = ["../features/pitch-yinjustin/"+mbid+".txt" for mbid in mbids if exists("../features/pitch-yinjustin/"+mbid+".txt")]

print filepaths
[an, abinCenters, acents] = iL.computeHist(filepaths)
apeakInfo = iL.findPeaks(an, abinCenters, lookahead=20, delta=0.00005, averageHist=True)

parameters = {}
for filepath in filepaths:
	[n, binCenters, cents] = iL.computeHist(filepath)
	peakInfo = iL.findPeaks(n, binCenters, lookahead=15, delta=0.00003, averageHist=False, refPeaks=apeakInfo["peaks"][0], refValleys=apeakInfo["valleys"][0])
	params = iL.characterizePeaks(peakInfo["peaks"][0], peakInfo["valleys"][0], cents)
	parameters[mbid] = params