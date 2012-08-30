# -*- coding: utf-8 -*-
"""
@author: Gopal
"""

import sys
import yaml
from os import listdir
from pprint import pprint
import DAO as dao

folder = sys.argv[1]
db = dao.DAO()
recordings = listdir(folder)
artists = {}
works = {}
vocal = 0
instru = 0
noInfo = 0

for i in recordings:
	mbid = i[:-4]
	recInfo = db.getRecordingInfo("mbid", mbid)
	if recInfo != "empty":
		if recInfo['artist']['name'] in artists.keys():
			artists[recInfo['artist']['name']].append(recInfo['title'])
		else:
			artists[recInfo['artist']['name']] = [recInfo['title']]
		vocal_flag = 0
		if 'artists' in recInfo.keys():
			for artist in recInfo['artists']:
				if artist['relation']['type'] == 'vocal':
					vocal_flag = 1
					vocal+=1
					break
			if not vocal_flag:
				instru+=1
		else:
			noInfo+=1
		if 'work' in recInfo.keys():
			if recInfo['work']['title'] in works.keys():
				works[recInfo['work']['title']].append(recInfo['artist']['name'])
			else:
				works[recInfo['work']['title']] = [recInfo['artist']['name']]

print "Folder: ", folder
#pprint(artists)
lens = [len(artists[i]) for i in artists.keys()]
if sum(lens) > len(lens):
	print "Per artist : ", lens
lens = [len(works[i]) for i in works.keys()]
if sum(lens) > len(lens):
	print "Per work : ", lens
#print 'Vocal & Instrumental:', vocal, "&", instru
print "----\n"
