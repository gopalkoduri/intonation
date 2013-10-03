# -*- coding:utf-8 -*-
#!/usr/bin/env python

from xml.dom import minidom
import xml.etree.ElementTree as ET
import numpy as np
import sys

cycleLength = 8

def redoAnnotationFile(annotationFile, resultFile):
    annotationCycles = []
    data = ET.parse(annotationFile)
    annotation = data.getroot()
    start = int(annotation[0][1][0].attrib['frame'])

    finalCycles = ET.Element('dataset', annotation[0][1].attrib)
    counter = 1
    
    #get the usual time-width between frames
    frameDiffs = []
    prevFrame = start
    for point in annotation[0][1][1:]:
        frameDiffs.append(int(point.attrib['frame'])-prevFrame)
        prevFrame = int(point.attrib['frame'])
    frameDiffs.sort()

    avgTime = np.mean(frameDiffs[:10])*1.2 #NOTE: Hardcoded

    for point in annotation[0][1][1:]:
        end = int(point.attrib['frame'])
        realEnd = end
        if (end-start) > avgTime:
            realEnd = end
            end = start+avgTime
        framePositions = np.linspace(start, end, cycleLength, endpoint=False)
        for i in xrange(len(framePositions)):
            newPoint =  ET.SubElement(finalCycles, 'point', {'label':\
                                        str(counter)+'.'+str(i%cycleLength+1), 'frame':str(int(framePositions[i]))})
        counter += 1
        start = realEnd #should be realEnd, and not just end.
    
    #the last annotation gets neglected!
    framePositions = np.linspace(start, start+avgTime, cycleLength, endpoint=False)
    for i in xrange(len(framePositions)):
        newPoint =  ET.SubElement(finalCycles, 'point', {'label':\
                                    str(counter)+'.'+str(i%cycleLength+1), 'frame':str(int(framePositions[i]))})

    annotation[0][1] = finalCycles
    res = ET.ElementTree(annotation)
    res.write(resultFile)

if __name__ == "__main__":
    redoAnnotationFile(sys.argv[1], sys.argv[2])
