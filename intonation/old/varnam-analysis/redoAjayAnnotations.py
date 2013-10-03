# -*- coding:utf-8 -*-
#!/usr/bin/env python

from xml.dom import minidom
import xml.etree.ElementTree as ET
import numpy as np
import sys

cycleLength = 8

def prettify(elem):
	"""Return a pretty-printed XML string for the Element.
	"""
	rough_string = ET.tostring(elem, 'utf-8')
	reparsed = minidom.parseString(rough_string)
	return reparsed.toprettyxml(indent="    ")

def redoAnnotationFile(annotationFile, resultFile):
	annotationCycles = []
	data = ET.parse(annotationFile)
	annotation = data.getroot()

	cycles = {}
	cycleSequence = []
	#annotation[0][1] is 'dataset' in annotations
	for point in annotation[0][1]:
		cycle = int(point.attrib['label'].split('.')[0])
		if cycle not in cycleSequence:
			cycles[cycle] = [point.attrib['frame']]
			cycleSequence.append(cycle)
		else:
			cycles[cycle].append(point.attrib['frame'])
	
	finalCycles = ET.Element('dataset', annotation[0][1].attrib)
	counter = 1
	for ind in cycleSequence:
		cycle = cycles[ind]
		if len(cycle) > cycleLength:
			numCycles = len(cycle)/cycleLength
			cycle = np.array(cycle)
			cycle = cycle.reshape(numCycles, cycleLength)
			for c in cycle:
				for i in xrange(len(c)):
					point = ET.SubElement(finalCycles, 'point', {'label': str(counter)+'.'+str(i+1), 'frame':str(c[i])})
				counter += 1
		else:
			for i in xrange(len(cycle)):
				point = ET.SubElement(finalCycles, 'point', {'label': str(counter)+'.'+str(i+1), 'frame':str(cycle[i])})
			counter += 1

	#return prettify(finalCycles)
	annotation[0][1] = finalCycles
	res = ET.ElementTree(annotation)
	res.write(resultFile)

if __name__ == "__main__":
	redoAnnotationFile(sys.argv[1], sys.argv[2])
