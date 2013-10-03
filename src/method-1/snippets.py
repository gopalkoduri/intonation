for rec in allParams.keys():
	for swara in allParams[rec].keys():
		for i in ["kurtosis", "skew2", "mean", "variance"]:
			allParams[rec][swara].pop(i)

def raaga(mbid):
	for r in raagaMBIDs.keys():
		if mbid in raagaMBIDs[r]:
			return r

import shutil
def concatFiles(files):
	dest = "-".join([i[:-5] for i in files])
	dest = dest+".arff"
	out = file("/homedtic/gkoduri/Workspace/intonationLib/data/method-2/experiments/experiment-alliedRaagas/experiment-allFeats/selectedSets/"+dest, "wb")
	shutil.copyfileobj(file("/homedtic/gkoduri/Workspace/intonationLib/data/method-2/experiments/data-allFeats/prelude-skew2.arff", "rb"), out)
	out.write("\n")
	out.write("@ATTRIBUTE segment {"+", ".join([i[:-5] for i in files])+"}\n\n@DATA\n")
	for f in files:
		shutil.copyfileobj(file("/homedtic/gkoduri/Workspace/intonationLib/data/method-2/experiments/data-allFeats/"+f, "rb"), out)
	out.close()

for f in randSets:
	files = f[:-5].split("-")
	files = [i+".arff" for i in files]
	concatFiles(files)

for group in ars:
	combos = combinations(group, 2)
	for files in combos:
		files = [f+".arff" for f in files]
		concatFiles(files)

for mbid in allParams.keys():
	r = raaga(mbid)
	if not r: continue
	r = r.replace(" ", "_")
	if not exists("data/"+r+"/"):
		mkdir("data/"+r+"/")
	yaml.dump(allParams[mbid], file("data/"+r+"/"+mbid+".sig", "w"), default_flow_style=False)


# filling up empty swaras and removing second higher octave
swaras = ["Sa_", "R1_", "R2/G1_", "G2/R3_", "G3_", "M1_", "M2_", "P_", "D1_", "D2/N1_", "D3/N2_", "N3_", "Sa", "R1", "R2/G1", "G2/R3", "G3", "M1", "M2", "P", "D1", "D2/N1", "D3/N2", "N3", "Sa^", "R1^", "R2/G1^", "G2/R3^", "G3^", "M1^", "M2^", "P^", "D1^", "D2/N1^", "D3/N2^", "N3^"]
others = ["Sa^^", "R1^^", "R2/G1^^", "G2/R3^^", "G3^^", "M1^^", "M2^^", "P^^", "D1^^", "D2/N1^^", "D3/N2^^", "N3^^"]

for rec in allParamsPhase2.keys():
    for swara in swaras:
        if swara not in allParamsPhase2[rec].keys():
            allParamsPhase2[rec][swara] = {"position":0, "mean":0, "variance":0, "skew1":0, "kurtosis":0, "amplitude":0, "skew2":0}
    for swara in allParamsPhase2[rec].keys():
        if swara in others:
            allParamsPhase2[rec].pop(swara)

savefig("bevahiour.pdf", dpi=300, facecolor='w', edgecolor='w', orientation='landscape', papertype=None, bbox_inches="tight", pad_inches=0.1)

count = 1
swaraDistributions = {}
for i in mbids:
	print count
	count += 1
	data = iL.pitch(i, vocal=True)
	sDs = iL.label_contours(data)
	swaraDistributions[i] = sDs
	
for f in listdir("."):
	data = file(f).readlines()
	data[0] = "@relation %s\n"%(f)
	out = file(f, "w")
	for line in data:
		out.write(line)
	out.close()


#dataset table
raagaMBIDs = yaml.load(file("raagaMBIDs-pruned.yaml"))
artists = {}
albums = {}
duration = {}
noexist = []
for raaga in raagaMBIDs.keys():
	for mbid in raagaMBIDs[raaga]:
		recInfo = db.getRecordingInfo("uuid", mbid)
		if recInfo == "empty":
			noexist.append(mbid)
			continue
		if raaga in albums.keys():
			albums[raaga].append(recInfo["release"]["title"])
		else:
			albums[raaga] = [recInfo["release"]["title"]]
		if raaga in artists.keys():
			artists[raaga].append(recInfo["artist"]["name"])
		else:
			artists[raaga] = [recInfo["artist"]["name"]]
		if raaga in duration.keys():
			duration[raaga] += recInfo["length"]
		else:
			duration[raaga] = recInfo["length"]

for raaga in raagas:
	print raaga.capitalize(), " & ", len(raagaMBIDs[raaga]), " & ", duration[raaga]/60000, " & ", len(unique(artists[raaga])), " & ", len(unique(albums[raaga])), "\\\\"
	

#feature selection statistics
paramStats = {}
counts = {}
dones = {}
count = len(stats)
for i in ["position", "amplitude", "mean", "variance", "skew2", "kurtosis"]:
	dones[i] = False
	counts[i] = 0

for i in stats.keys():
	for d in dones.keys():
		dones[d] = False
	for param in stats[i].keys():
		if not dones[param]:
			counts[param]+=1
			dones[param] = True
		if param in paramStats.keys():
			paramStats[param] += stats[i][param]
		else:
			paramStats[param] = stats[i][param]
			

#stability plots
data = loadtxt("results-clean.csv", delimiter=",", dtype="float", usecols=[2,4,6,8,10,12])
classifiers = ["naive bayes", "k-nearest neighbours", "support vector machine", "random forest", "logistic regression", "multilayer peceptron"]

linestyles = ["-.", "-", "--"]
markers = ["*", "x", "d"]

for i in xrange(len(data[0])):
	count = 1
	_sum = 0
	x = []
	y = []
	for j in data[:, i]:
		x.append(count)
		_sum+=j
		y.append(1.0*_sum/count)
		count+=1
	x = x[::20]
	y = y[::20]
	if i < len(linestyles):
		plot(x, y, lw=2.0, linestyle=linestyles[i], color="k")
	else:
		plot(x, y, lw=2.0, linestyle=':', marker=markers[i-len(linestyles)], color="k")
legend(classifiers, prop={"size":18})
xticks(fontsize=20)
yticks(fontsize=20)
xlabel("Number of experiments", fontsize=20)
ylabel("Accuracy", fontsize=20)
