import intonationLib as iL
import numpy as np
import pickle
import yaml
from scipy.stats import skew, variation, kurtosis, mode
from os.path import exists
from os import mkdir

#swaras = ["Sa_", "R1_", "R2/G1_", "G2/R3_", "G3_", "M1_", "M2_", "P_", "D1_", "D2/N1_", "D3/N2_", "N3_", "Sa", "R1", "R2/G1", "G2/R3", "G3", "M1", "M2", "P", "D1", "D2/N1", "D3/N2", "N3", "Sa^", "R1^", "R2/G1^", "G2/R3^", "G3^", "M1^", "M2^", "P^", "D1^", "D2/N1^", "D3/N2^", "N3^", "Sa^^", "R1^^", "R2/G1^^", "G2/R3^^", "G3^^", "M1^^", "M2^^", "P^^", "D1^^", "D2/N1^^", "D3/N2^^", "N3^^"]
#swaras = ["Sa", "R1", "R2/G1", "G2/R3", "G3", "M1", "M2", "P", "D1", "D2/N1", "D3/N2", "N3", "Sa^", "R1^", "R2/G1^", "G2/R3^", "G3^", "M1^", "M2^", "P^", "D1^", "D2/N1^", "D3/N2^", "N3^"]
swaras = ["D1", "D1^"]

allSDs = yaml.load(file("../../data/allParameters-25-50.yaml", "r"))
raagas = ['begada', 'bhairavi', 'hindolam', 'kambhoji', 'manji', 'mukhari', 'thodi']
#raagas = ['bhairavi', 'hindolam', 'thodi']
allParams = {}

for raaga in raagas:
	allParams[raaga] = {}
	temp = yaml.load(file("../../data/exp1/"+raaga+".yaml", "r"))
	mbids = temp.keys()
	for mbid in mbids:
		if mbid not in allParams[raaga].keys(): continue
		parameters = {}
		_total = sum([len(allSDs[mbid][swara]) for swara in allSDs[mbid].keys()])
		for swara in swaras:
			if swara in allSDs[mbid].keys() and len(allSDs[mbid][swara]) > 0:
				_mode = mode(allSDs[mbid][swara].astype(int))
				p = float(_mode[0][0])
				amplitude = float(_mode[1][0]*1.0/_total)
				_mean = float(np.mean(allSDs[mbid][swara]))
				_variance = float(variation(allSDs[mbid][swara]))
				_skew = float(skew(allSDs[mbid][swara]))
				_kurtosis = float(kurtosis(allSDs[mbid][swara]))
				if _variance != 0:
					pearsonSkew = float(3.0*(_mean-p)/np.sqrt(abs(_variance)))
				else:
					pearsonSkew = 0
				parameters[swara] = {"position":float(p), "mean":_mean, "variance":_variance, "skew1":_skew, "kurtosis":_kurtosis, "amplitude":amplitude, "skew2":pearsonSkew}
			else:
				parameters[swara] = {"position":0, "mean":0, "variance":0, "skew1":0, "kurtosis":0, "amplitude":0, "skew2":0}

		#allParams[raaga][mbid] = parameters
		if not exists("../../data/exp1/weka/"+raaga+"/"):
			mkdir("../../data/exp1/weka/"+raaga)
		else:
			yaml.dump(parameters, file("../../data/exp1/weka/"+raaga+"/"+mbid+".sig", "w"))

#yaml.dump(allParams, file("../../data/exp1/allParameters.yaml", "w"), default_flow_style=False)
