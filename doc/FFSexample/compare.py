import sys

sys.path.append('../../libs')

import argparse
import time
from math import *

sys.path.append('../../')

from metaclass import twiss
from mapclass import Map2

### Argument reading

parser = argparse.ArgumentParser(description="FFS test for MapClass2")
parser.add_argument('-o', help="Order at which run the test. Default is 3 (for fast results), 6 would be in the order of 3 minutes.",
                    type=int, default=3, dest='order')
parser.add_argument('-g', help="Run with Gaussian Delta or not. By default it doesn't",
                    action='store_true', dest='gaussian')

args = parser.parse_args()

### Setup

startTime = time.time()

t = twiss('assets/ffs.twiss')

m = Map2(t,order=args.order)

mm = Map2(filename='assets/fort.18',order=args.order)

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
