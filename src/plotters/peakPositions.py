import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Polygon
import yaml

allParameters = yaml.load(file("../../data/allParameters-25-50.yaml"))

param = "position"
paramValues = {}
#swaras = ["R1", "R2/G1", "G2/R3", "G3", "M1", "M2", "P", "D1", "D2/N1", "D3/N2", "N3", "Sa^"]
#swaras = ["Sa", "R1", "R2/G1", "G2/R3", "G3", "M1", "M2", "P", "D1", "D2/N1", "D3/N2", "N3", "Sa^", "R1^", "R2/G1^", "G2/R3^", "G3^", "M1^", "M2^", "P^", "D1^", "D2/N1^", "D3/N2^", "N3^"]

#raagas = ["bhairavi", "saveri", "anandabhairavi", "khamas", "thodi",  "pantuvarali", "kalyani", "sourashtram"]
#swaras = ["P"]
#templatePeak = 700.18

#raagas = ["anandabhairavi", "khamas", "hindolam", "sourashtram"]
#raagas.sort(reverse=False)
#swaras = ["D3/N2"]
#templatePeak = 987.87

#raagas = ["khamas", "kalyani", "sourashtram"]
#swaras = ["D2/N1"]
#templatePeak = 889.58

raagas = ["bhairavi", "saveri", "anandabhairavi", "khamas", "thodi",  "hindolam", "sourashtram"]
swaras = ["M1"]
templatePeak = 495.05

#raagas = ["saveri", "khamas", "pantuvarali", "kalyani"]
#swaras = ["G3"]
#raagas = ["bhairavi", "anandabhairavi", "kalyani", "mukhari", "kambhoji"]
#swaras = ["R2/G1"]


foundSwaras = []

for raaga in raagas:#allParameters.keys():
	paramValues[raaga] = {}
	for rec in allParameters[raaga].keys():
		for swara in allParameters[raaga][rec].keys():
			if swara in swaras:
				if swara not in foundSwaras: foundSwaras.append(swara)
				if swara not in paramValues[raaga].keys():
					paramValues[raaga][swara] = [allParameters[raaga][rec][swara][param]]
				else:
					paramValues[raaga][swara].append(allParameters[raaga][rec][swara][param])

data = []

#removeSwaras = ["Sa"]
#for i in removeSwaras:
	#foundSwaras.remove(i)

foundSwaras.sort()
for swara in foundSwaras:
	for i in xrange(len(raagas)):
		if swara not in paramValues[raagas[i]].keys():
			data.append([0])
		else:
			data.append(paramValues[raagas[i]][swara])
			#print raagas[i], swara, paramValues[raagas[i]][swara]

#numSwaras = len(paramValues[raagas[0]].keys())
numSwaras = len(foundSwaras)

axisFontSize=40
swaraLabelFontSize=28
medianValueFontSize=28

fig = plt.figure(figsize=(10,6))
fig.canvas.set_window_title(param+" of "+",".join(raagas))
ax1 = fig.add_subplot(111)
plt.subplots_adjust(left=0.095, right=0.95, top=0.9, bottom=0.25)

bp = plt.boxplot(data, notch=0, sym='+', vert=1, whis=1.5)
bottom, top = plt.ylim()
plt.setp(bp['boxes'], color='black')
plt.setp(bp['whiskers'], color='black')
plt.setp(bp['fliers'], color='black', marker='+')
plt.axhline(templatePeak, color="black", ls="--", lw=2.5)

# Add a horizontal grid to the plot, but make it very light in color
# so we can use it for reading data values but not be distracting
ax1.yaxis.grid(True, linestyle='-', which='major', color='lightgrey',
              alpha=0.5)

# Hide these grid behind plot objects
ax1.set_axisbelow(True)
_title1 = param.title()+" of "
_title2 = ", ".join(raagas)
#ax1.set_title(_title1+_title2.title())
#ax1.set_xlabel("M1", fontsize=axisFontSize)
ax1.set_ylabel('Cents', fontsize=axisFontSize)

# Now fill the boxes with desired colors
#boxColors = ['darkkhaki',"#9F8170", "#C19A6B"]
#boxColors = ["#3C1414", "#960018", "#00416A", "#87A96B", "#CD9575", "#873260", "#79443B", "#B5A642", "#0093AF", "#004225"]
boxColors = np.repeat("#777777", len(raagas))
#boxColors = boxColors[:len(raagas)]
numBoxes = numSwaras*len(raagas)
medians = range(numBoxes)
for i in range(numBoxes):
	box = bp['boxes'][i]
	boxX = []
	boxY = []
	for j in range(5):
		boxX.append(box.get_xdata()[j])
		boxY.append(box.get_ydata()[j])
	boxCoords = zip(boxX,boxY)
	# Alternate between colors
	k = i % len(raagas)
	boxPolygon = Polygon(boxCoords, facecolor=boxColors[k])
	ax1.add_patch(boxPolygon)
	# Now draw the median lines back over what we just filled in
	med = bp['medians'][i]
	medianX = []
	medianY = []
	for j in range(2):
		medianX.append(med.get_xdata()[j])
		medianY.append(med.get_ydata()[j])
		plt.plot(medianX, medianY, 'k')
		medians[i] = medianY[0]
  # Finally, overplot the sample averages, with horixzontal alignment
  # in the center of each box
	plt.plot([np.average(med.get_xdata())], [np.average(data[i])],
           color='w', marker='*', markeredgecolor='k')

# Set the axes ranges and axes labels
#ax1.set_xlim(0.5, numBoxes+0.5)
#top = top*1.1
#bottom = bottom*1.1
#ax1.set_ylim(bottom, top)
xtickNames = plt.setp(ax1, xticklabels=raagas)
plt.setp(xtickNames, rotation=30, fontsize=swaraLabelFontSize)

# Due to the Y-axis scale being different across samples, it can be
# hard to compare differences in medians across the samples. Add upper
# X-axis tick labels with the sample medians to aid in comparison
# (just use two decimal places of precision)

pos = np.arange(numBoxes)+1
if max(medians) > 0.01:
	upperLabels = [str(np.round(s, 2)) for s in medians]
else:
	upperLabels = [str(np.round(s, 2)) for s in medians]
#weights = ['bold', 'semibold', 'normal', 'bold', 'semibold', 'normal', 'bold', 'semibold', 'normal', 'bold', 'semibold', 'normal']
weights = np.repeat("normal", len(raagas))
for tick,label in zip(range(numBoxes),ax1.get_xticklabels()):
	k = tick % len(raagas)
	ax1.text(pos[tick], top*.996, upperLabels[tick],
			horizontalalignment='center', size=medianValueFontSize, weight=weights[k],
			color=boxColors[k])

#Finally, add a basic legend
x = 0.8
y = 0.8
ystep = 0.04
xstep = 0.03

#plt.figtext(x, y, '--', color='yellow', backgroundcolor='#960018', weight='roman', size='26')
#plt.figtext(x, y-0.03, '  Mean from general template', color='black', weight='roman', size=swaraLabelFontSize)
plt.show()