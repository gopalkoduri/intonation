# -*- coding:utf-8 -*-
#/usr/bin/env python

import sys
import yaml
from os import system

def chunks(l, n):
    """ Yield successive n-sized chunks from l.
    """
    for i in xrange(0, len(l), n):
        yield l[i:i+n]

label = sys.argv[1]
fingerLen = 4

data = file(sys.argv[2], 'r').readlines()
notation = []
for line in data:
	cycle = line.split()
	cycle = [i.strip('-') for i in cycle]
	cycle = list(chunks(cycle, fingerLen))
	notation.append(cycle)
	
varnam  = {label:notation}
yaml.dump(varnam, file(sys.argv[3], 'w'), default_flow_style=False)
#system('cat %s'%(sys.argv[3]))
