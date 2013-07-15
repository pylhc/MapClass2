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
parser.add_argument('-o', help="Order at which run the test. Default is 3",
                    type=int, default=3, dest='order')
parser.add_argument('-g', help="Run with Gaussian Delta or not. By default it doesn't",
                    action='store_true', dest='gaussian')
parser.add_argument('-n', help="How many cores to run on. By default it is 1",
                    type=int, default=1, dest='nbProc')                    
parser.add_argument('-m', help="Which formula to use for multipoles",
                    type=int, default=0, dest='fmultipole')
parser.add_argument('-ns', help="Do not strip the line of unused elements (MATRIX, MONITOR, MARKER): default true",
                    action='store_false', default=True, dest='strpl')
                    

args = parser.parse_args()

### Setup

startTime = time.time()

t = twiss2('assets/ffs.twiss')
#t = t.stripLine(); already done in c++ by default
# The following is a call to construct the Map using the C++ Map construction from a twiss file
#m = Map2('assets/ffs.twiss',filenameerr=None, order=args.order, nbProc=args.nbProc, strpl=args.strpl) 

# The following is call to the old mapclass which uses Python map construction
#m = Map2(t, "old", order=args.order) 

# The following is a call to construct the Map from a Twiss Python object
m = Map2(t, terr=None, order=args.order, nbProc=args.nbProc, fmultipole=args.fmultipole, strpl=args.strpl)

mm = Map2(filename='assets/fort.18',order=args.order) #from fort call

v = ['x','px','y','py']

### Obtain numbers

print ''
print "!!!!! Maps of order %i" % args.order,

if args.gaussian:
    print "and with Gaussian Delta"
print ''
print "########## Comparison [It's expected of both to give the same value for high chi2]"

print "* MapClass.comp() #It ignores coefficients where 'd' isn't 0." 
print "=> %e" % mm.comp(m,v).real
print ''
print "* MapClass.compc() #It compares all the elements."
print "=> %e" % mm.compc(m, v).real

print ''
print "[T] Represents the value calculated from the Twiss where we do the calculations to get the Map."
print "[F] Represents the value calculated from the coefficients obtained by PTC."
print ''

###                                                                                                     
betx=64.99988501
bety=17.99971417
gamma=3e6
ex=68e-8
ey=2e-8
sigmaFFS=[sqrt(ex*betx/gamma), sqrt(ex/betx/gamma), sqrt(ey*bety/gamma), sqrt(ey/bety/gamma), 0.01, 0.002]
###

print "@@@@@@@@ Sigma"

for i in v:
    print "[T][%s] %e" % (i, sqrt(m.sigma(i,sigmaFFS,args.gaussian).real))
    print "[F][%s] %e" % (i, sqrt(mm.sigma(i,sigmaFFS,args.gaussian).real))

print ''
print "@@@@@@@@ Offset"

for i in v:
    print "[T][%s] %e" % (i, m.offset(i,sigmaFFS).real)
    print "[F][%s] %e" % (i, mm.offset(i,sigmaFFS).real)


print ''
print 'It took', time.time() - startTime, 'seconds to run.'   
