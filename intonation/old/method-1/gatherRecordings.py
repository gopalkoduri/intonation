#!/usr/bin/env python

import sys
sys.path.append("../")
import yaml

import DAO as dao
import intonationLib as iL

db = dao.DAO()
artists = db.getArtists(collection='carnatic')

recData = {}
for artist in artists[1]:
	#Get all recordings sorted by raaga, for each artist.
	artistRecs = db.getRecordingsByArtist('mbid', artist['uuid'])
	artistRecsByRaaga = {}
	for rec in artistRecs[1]:
		recInfo = db.getRecordingInfo('mbid', rec['uuid'])
		if 'tags' in recInfo.keys():
			for i in recInfo['tags']:
				if i['category'] == 'raaga':
					if i['tag'] in artistRecsByRaaga.keys():
						artistRecsByRaaga[i['tag']].append({'uuid':recInfo['uuid'], 'title':recInfo['title'], 'raaga':i['tag']})
					else:
						artistRecsByRaaga[i['tag']] = [{'uuid':recInfo['uuid'], 'title':recInfo['title'], 'raaga':i['tag']}]
	
	#Handle spelling mistakes of raagas
	temp = {}
	raagas = artistRecsByRaaga.keys()
	print raagas
	analyzed = []
	for raaga in raagas:
		temp[raaga] = artistRecsByRaaga[raaga]
		duplicates = iL.raagaDuplicates(raaga)
		for duplicate in duplicates:
			if raaga == duplicate: continue
			if (duplicate not in analyzed) and (duplicate in raagas):
				analyzed.append(duplicate)
				temp[raaga].extend(artistRecsByRaaga[duplicate])
	artistRecsByRaaga = temp
	
	#Keep only those raagas which have at least n recordings per artist.
	n = 2
	for raaga in artistRecsByRaaga.keys():
		if len(artistRecsByRaaga[raaga]) < n:
			artistRecsByRaaga.pop(raaga)
	recData[artist['uuid']] = {'name':artist['name'], 'recordings':artistRecsByRaaga}

#clear recData of artists who have empty recordings
for artistID in recData.keys():
	if len(recData[artistID]['recordings']) == 0:
		recData.pop(artistID)

#Finally retain only those raagas which have at least n artists
n = 2
retainedRaagas = {}
for artistID in recData.keys():
	for raaga in recData[artistID]['recordings'].keys():
		if raaga in retainedRaagas.keys():
			retainedRaagas[raaga] += 1
		else:
			retainedRaagas[raaga] = 1

for raaga in retainedRaagas.keys():
	if retainedRaagas[raaga] < n:
		retainedRaagas.pop(raaga)

for artistID in recData.keys():
	for raaga in recData[artistID]['recordings'].keys():
		if raaga not in retainedRaagas.keys():
			recData[artistID]['recordings'].pop(raaga)

#clear recData of artists who have empty recordings
vocalists = yaml.load(file("vocalists.yaml"))
for artistID in recData.keys():
	if (len(recData[artistID]['recordings']) == 0) or (artistID not in vocalists.keys()):
		recData.pop(artistID)

yaml.dump(recData, file('vocal-raaga-artist-constrained.yaml', 'w'))