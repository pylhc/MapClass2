#! /usr/bin/env python

import sys
import time
from math import *

if len(sys.argv) < 3:
  sys.stderr.write('Usage: sys.argv[0] dir [-t|-f]')
  sys.exit(1)

sys.path.append(sys.argv[1])

from metaclass2 import twiss2
from mapclass import Map2

o = int(sys.argv[3])

startTime = time.time()

if sys.argv[2] == "-t":
  t = twiss2()
  m = Map2(t,order=o)
if sys.argv[2] == "-m":
  t = twiss2()
  t = t.stripLine()
  t = t.mergeElems()
  m = Map2(t,order=o)
else:
  m = Map2(order=o)

v = ['x','px','y','py']

###
betx=64.99988501
bety=17.99971417
gamma=3e6
ex=68e-8
ey=2e-8
sigmaFFS=[sqrt(ex*betx/gamma), sqrt(ex/betx/gamma), sqrt(ey*bety/gamma), sqrt(ey/bety/gamma), 0.01, 0.002]
###

for i in v:
    m.sigma(i,sigmaFFS)

print time.time() - startTime
