import yaml
import sys
sys.path.append("../")
import DAO as dao
from shutil import copy

db = dao.DAO()
recData = yaml.load(file("vocal-raaga-artist-constrained.yaml"))
pathData = yaml.load(file("/media/Data/Music/CompMusic/Carnatic/metadata/Carnatic.yaml"))
notFound = []

for artistID in recData.keys():
	for raaga in recData[artistID]['recordings'].keys():
		for rec in recData[artistID]['recordings'][raaga]:
			if rec['uuid'] in pathData.keys():
				#copy(pathData[rec['uuid']]['path'], 'data/'+rec['uuid']+".mp3")
				print pathData[rec['uuid']]['path']
			else:
				notFound.append(rec['uuid'])

#print "Paths not found for ", len(notFound), "songs."
#for i in notFound:
#	recInfo = db.getRecordingInfo('mbid', i)
#	print i, recInfo['release']['title']
