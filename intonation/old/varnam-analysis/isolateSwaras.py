# -*- coding:utf-8 -*-
#!/usr/bin/env python

import yaml
import xml.etree.ElementTree as ET
import intonationLib as iL
import numpy as np


cycleLength = 8
frameRate = 44100.0

def isolateSwaras(notationFile, pitchFile, annotationFile):
    notation = yaml.load(file(notationFile))
    pitch = np.loadtxt(pitchFile, delimiter='\t', usecols=[0, 1])

    annotationCycles = []
    data = ET.parse(annotationFile)
    annotation = data.getroot()
    prevCycle = 1
    start = float(annotation[0][1][0].attrib['frame'])
    prevFrame = start

    #annotation[0][1] is 'dataset' in annotations
    for point in annotation[0][1]:
        cycle = int(point.attrib['label'].split('.')[0])
        if cycle > prevCycle:
            #divide cycleLength-1 since prevFrame corresponds to n.7 in adi taala
            frameWidth = (prevFrame-start)/(cycleLength-1)
            end = prevFrame+frameWidth
            annotationCycles.append([start/frameRate, end/frameRate])
            start = float(point.attrib['frame'])
            prevCycle = cycle
        else:
            prevFrame = float(point.attrib['frame'])
    #the last cycle gets neglected.. 
    frameWidth = (prevFrame-start)/(cycleLength-1)
    end = prevFrame+frameWidth
    annotationCycles.append([start/frameRate, end/frameRate])
    
    #Pseudo code
    #
    #get the sequence of 'cycles of notation' based on the length of 
    #each section and its repetition (some sections are sung once, some twice)
    #
    #for each cycle in the sequence:
    #    get timestamps of start and end of the cycle from annotation
    #    get pitch of that timestamp and cut it into 32/64 pieces
    #    push each piece into the bin of corresponding svara
    
    fastCycles = []
    notationCycles = []
    notationCycles.extend(notation['pallavi'])
    notationCycles.extend(notation['anupallavi'])
    notationCycles.extend(notation['muktayiswaram'])

    fastCycles.extend(len(notationCycles)+np.array(range(len(notation['pallavi'])+\
            len(notation['anupallavi'])+len(notation['muktayiswaram']))))
    notationCycles.extend(notation['pallavi'])
    notationCycles.extend(notation['anupallavi'])
    notationCycles.extend(notation['muktayiswaram'])
    
    #print fastCycles

    for cswaram in notation['chittiswaram']:
        #first speed
        notationCycles.extend(notation['charanam'])
        notationCycles.extend(cswaram)
        #second speed
#        fastCycles.extend(len(notationCycles)+np.array(range(len(notation['charanam'])+len(cswaram))))
#        notationCycles.extend(notation['charanam'])
        if len(cswaram) % 2 == 1:
            #print "counting once"
            fastCycles.extend(len(notationCycles)+np.array(range(len(notation['charanam'])+len(cswaram))))
            notationCycles.extend(notation['charanam'])
        else:
            #print "counting twice"
            fastCycles.extend(len(notationCycles)+np.array(range(2*len(notation['charanam'])+len(cswaram))))
            notationCycles.extend(notation['charanam'])
            notationCycles.extend(notation['charanam'])
        notationCycles.extend(cswaram)
    #And then it ends with one last cycle of charanam
    notationCycles.extend(notation['charanam'])
    #return notationCycles
    #return fastCycles

    distributions = {}
    placeCount = 0
    startFlag = 0


    #if slow cycle:
    #    2 cycles in annotation = 1 cycle in notation
    #if fast cycle:
    #    1 cycle in annotation = 1 cycle in notation
    #annInd and ind progress by these rules

    annInd = 0
    ind = 0
    annLayerPoints = ""
    annLayerStart = -1

    prevSwara = notationCycles[0][0][0]
    while ind < len(notationCycles):
        swaras = np.array(notationCycles[ind])
        numSwaras = len(notationCycles[ind])*len(notationCycles[ind][0])
        swaras = list(swaras.reshape(numSwaras))
        #replace commas with previous swara label
        for i in xrange(len(swaras)):
            if swaras[i] == '': print notationCycles[ind]
            if swaras[i] == ',':
                swaras[i] = prevSwara
            else:
                prevSwara = swaras[i]

        #If it is a slow cycle, add another cycle from notation so that the whole thing becomes one cycle in annotation
        annStart = annotationCycles[annInd][0]
        annEnd = annotationCycles[annInd][1]
        if ind not in fastCycles:
            annInd = annInd+1
            annEnd = annotationCycles[annInd][1]
        if annLayerStart == -1:
            annLayerStart = annStart

        #NOTE: Make sure the cycles in annotation which does not correspond to notation are manually removed
        annotationDivisions = np.linspace(annStart, annEnd, numSwaras+1)
        for i in xrange(numSwaras):
            startIndex = iL.findNearestIndex(pitch[:, 0], annotationDivisions[i])
            endIndex = iL.findNearestIndex(pitch[:, 0], annotationDivisions[i+1])
            annLayerPoints = annLayerPoints+\
                    """<point frame="%d" height="0.6" label="%s"/>\n"""\
                    %(int(frameRate*pitch[startIndex, 0]), swaras[i])
            #annLayer = annLayer+\
            #        """%f 0.8   %s\n"""%(pitch[startIndex, 0], swaras[i])
            #print annotationDivisions[i], annotationDivisions[i+1], startIndex, endIndex
            if swaras[i] in distributions.keys():
                distributions[swaras[i]].extend(pitch[:, 1][startIndex:endIndex])
            else:
                distributions[swaras[i]] = list(pitch[:, 1][startIndex:endIndex])

        ind = ind+1
        annInd+=1
    
    annLayerEnd = annEnd
    annLayerPrelude = """<?xml version="1.0" encoding="UTF-8"?>
    <!DOCTYPE sonic-visualiser>
    <sv>
      <data>
          <model id="1" name="" sampleRate="44100" start="%d" end="%d"
          type="sparse" dimensions="2" resolution="1" notifyOnAdd="true"
          dataset="0"  subtype="text"/>
              <dataset id="0" dimensions="2">

    """%(annLayerStart, annLayerEnd)
    annLayerConclusion = """
                </dataset>
          </data>
            <display>
                <layer id="4" type="text" name="Text &lt;2&gt;" model="1"
                colourName="Red" colour="#800000" darkBackground="false" />
          </display>
      </sv>
    """
    annLayer = annLayerPrelude+annLayerPoints+annLayerConclusion

    print "notation cycles:", ind,"/",len(notationCycles)
    print "annotation cycles:", annInd,"/",len(annotationCycles)
    for swara in distributions.keys():
        distributions[swara] = np.array(distributions[swara])
        distributions[swara] = distributions[swara][distributions[swara] != -np.inf]
    return [annLayer, distributions]
