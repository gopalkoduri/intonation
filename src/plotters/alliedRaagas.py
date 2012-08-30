import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Polygon
import yaml

allParameters = yaml.load(file("../../data/allParameters-25-50.yaml"))

param = "skew2"
_ylabel = "Skewness"
paramValues = {}
#swaras = ["R1", "R2/G1", "G2/R3", "G3", "M1", "M2", "P", "D1", "D2/N1", "D3/N2", "N3", "Sa^"]
#swaras = ["Sa", "R1", "R2/G1", "G2/R3", "G3", "M1", "M2", "P", "D1", "D2/N1", "D3/N2", "N3", "Sa^", "R1^", "R2/G1^", "G2/R3^", "G3^", "M1^", "M2^", "P^", "D1^", "D2/N1^", "D3/N2^", "N3^"]

raagas = ["mukhari", "bhairavi", "manji"]
#Kurtosis - mukhari, bhairavi, manji
#swaras = ["R2/G1", "M1", "Sa^", "R2/G1^"]
#Skew - mukhari, bhairavi, manji
swaras = ["D2/N1", "R2/G1", "P", "P^"]

#raagas = ["begada", "kambhoji"]
#Variance - Begada, Kambhoji
#swaras = ["G3", "G3^", "R2/G1", "R2/G1^"]
#Kurtosis - Begada, Kambhoji
#swaras = ["G3", "G3^", "D2/N1", "D2/N1^"]

#raagas = ["suratti", "kedaragowla"]
##Skew
#swaras = ["M1", "M1^", "P", "P^"]

#raagas = ["bhairavi", "saveri", "anandabhairavi", "khamas", "thodi",  "pantuvarali", "kalyani", "sourashtram"]
#swaras = ["P"]

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

foundSwaras.sort(reverse=False)
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

# Add a horizontal grid to the plot, but make it very light in color
# so we can use it for reading data values but not be distracting
ax1.yaxis.grid(True, linestyle='-', which='major', color='lightgrey',
              alpha=0.5)

# Hide these grid behind plot objects
ax1.set_axisbelow(True)
_title1 = param.title()+" of "
_title2 = ", ".join(raagas)
#ax1.set_title(_title1+_title2.title())
#ax1.set_xlabel('Swaras', fontsize=axisFontSize)
ax1.set_ylabel(_ylabel, fontsize=axisFontSize)

# Now fill the boxes with desired colors
#boxColors = ['darkkhaki',"#9F8170", "#C19A6B"]
#boxColors = ["#000000", "#cccccc", "#333333", "#eeeeee", "#666666", "#999999", "#79443B", "#B5A642", "#0093AF", "#004225"]
#boxColors = ["#000000", "#999999"]
boxColors = ["#000000", "#777777", "#cccccc"]
#boxColors = ["#000000", "#555555", "#999999", "#cccccc"]

boxColors = boxColors[:len(raagas)]
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
           color='w', marker='*', markeredgecolor='k', label=i)

# Set the axes ranges and axes labels
#ax1.set_xlim(0.5, numBoxes+0.5)
ax1.set_ylim(bottom, top*1.1)
xtickNames = plt.setp(ax1, xticklabels=np.repeat(foundSwaras, len(raagas)))
plt.setp(xtickNames, rotation=45, fontsize=swaraLabelFontSize)

# Due to the Y-axis scale being different across samples, it can be
# hard to compare differences in medians across the samples. Add upper
# X-axis tick labels with the sample medians to aid in comparison
# (just use two decimal places of precision)

pos = np.arange(numBoxes)+1
upperLabels = [str(np.round(s, 1)) for s in medians]
weights = ['bold', 'semibold', 'normal', 'bold', 'semibold', 'normal', 'bold', 'semibold', 'normal', 'bold', 'semibold', 'normal']
for tick,label in zip(range(numBoxes),ax1.get_xticklabels()):
	k = tick % len(raagas)
	ax1.text(pos[tick], top, upperLabels[tick],
			horizontalalignment='center', size=medianValueFontSize, weight=weights[k],
			color=boxColors[k])

#Finally, add a basic legend
x = 0.75
y = 0.4
ystep = 0.04
xstep = 0.03
for i in xrange(len(raagas)):
	plt.figtext(x, y+0.02,  "__", backgroundcolor=boxColors[i], color=boxColors[i], weight='roman', size=swaraLabelFontSize)
	plt.figtext(x+xstep, y,  raagas[i], color='black', weight='roman', size=swaraLabelFontSize)
	y = y-ystep

plt.figtext(x, y, '*', color='white', backgroundcolor='black', weight='roman', size='medium')
plt.figtext(x+xstep, y, 'mean', color='black', weight='roman', size=swaraLabelFontSize)

plt.show()