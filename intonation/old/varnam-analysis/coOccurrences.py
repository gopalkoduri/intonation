#/usr/bin/env python

import yaml
import pickle
import numpy as np

def swaraSequence(notationFile):
    notation = yaml.load(file(notationFile))
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
    return np.concatenate(np.concatenate(notationCycles))


if __name__ == "__main__":
    raagas = ["abhogi", "begada", "mohanam", "kalyani", "sahana", "saveri", "shree"]
    for r in raagas:
        swaras = swaraSequence(r+".yaml")

        prevSwara = swaras[0]
        i = 1
        while i < len(swaras):
            if swaras[i] == '': print r, i
            if swaras[i] == ',':
                swaras[i] = prevSwara
            else:
                prevSwara = swaras[i]
            i += 1

        uniqSwaras = list(np.unique(swaras))
        countMatrix = [[0 for x in xrange(len(uniqSwaras))] for x in xrange(len(uniqSwaras))]

        prevSwara = swaras[0]
        i = uniqSwaras.index(prevSwara)
        for swara in swaras[1:]:
            j = uniqSwaras.index(swara)
            countMatrix[i][j] += 1
            i = j
        pickle.dump([uniqSwaras, countMatrix], file(r+"-countMatrix.pickle", "w"))
