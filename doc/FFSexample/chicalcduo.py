import sys
import cProfile
import pstats

from ctypes import *
cdll.LoadLibrary("../../libs/boost_1_53_0/libboost_python.so.1.53.0")

sys.path.append('../../libs')

import argparse
import time
from math import *

sys.path.append('../../')

from metaclass2 import twiss2
from mapclass import Map2
#from mapclassparallel import Map2
#from mapclasspool import Map2


### Argument reading

parser = argparse.ArgumentParser(description="FFS test for MapClass2")
parser.add_argument('-o', help="Order at which run the test. Default is 3 (for fast results), 6 would be in the order of 3 minutes.",
                    type=int, default=3, dest='order')
parser.add_argument('-g', help="Run with Gaussian Delta or not. By default it doesn't",
                    action='store_true', dest='gaussian')
parser.add_argument('-n', help="How many cores to run on. By default it is 1",
                    type=int, default=1, dest='nbProc')                    
                    

args = parser.parse_args()

t = twiss2('assets/ffs.twiss')
t = t.stripLine()
t = t.mergeElems()

sext1 = [24, 38, 43, 74, 82]
sext2 = [24, 30, 31, 32, 38, 43, 49, 74, 82]

#for i in range(0, len(sext1)):
#    print t.elems[sext1[i]].S
    
positions = [0.07967687682559, 0.049979477525608, 0.052613241375277, 0.193170454943335, 0.00738902713769]
strengths = [10.2881885328145, 133.508521888545, -109.469292468306, 68.8207501444652, -5.09437310259855, 15.5187207629783, -1.27849821750552, -6.01584279027468, 21.8375550635217]

#sys.stdout.write(".")
sys.stdout.flush()

for i in range(0, len(positions)):

    t0 = t.alterElem(sext1[i], dPos = positions[i])
    if t is None: 
        print "IT NO WORK"
        break
    else: 
        t = t0
#        print t.elems[sext1[i]].S
    
for i in range(0, len(strengths)):
    t.elems[sext2[i]].K2L = strengths[i]
    #print t.elems[sext2[i]].K2L
         
m = Map2(t, terr=None, order=args.order, nbProc=args.nbProc)
###                                                                                                     
betx=64.99988501
bety=17.99971417
gamma=3e6
ex=68e-8
ey=2e-8
sigmaFFS=[sqrt(ex*betx/gamma), sqrt(ex/betx/gamma), sqrt(ey*bety/gamma), sqrt(ey/bety/gamma), 0.01, 0.002]
###
h = sqrt(m.sigma('x', sigmaFFS, args.gaussian).real)
j = sqrt(m.sigma('y', sigmaFFS, args.gaussian).real)
h0 = 40e-9 #sigmaFFS[0]
j0 = 1e-9 #sigmaFFS[2]

chi = (h-h0)**2/40**2 + (j-j0)**2
print 'chi2 = ',chi, 'sigx = ', h, 'sigy = ', j

