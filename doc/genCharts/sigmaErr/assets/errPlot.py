import sys

from math import *

if len(sys.argv) < 2:
  sys.stderr.write('Usage: sys.argv[0] dir')
  sys.exit(1)

path = sys.argv[1]

sys.path.append(path)

from metaclass2 import twiss2
from mapclass import Map2

### Setup

class DummyFile(object):
    def write(self, x): pass

###                                                                                                     
v = ['x','px','y','py']

betx=64.99988501
bety=17.99971417
gamma=3e6
ex=68e-8
ey=2e-8
sigmaFFS=[sqrt(ex*betx/gamma), sqrt(ex/betx/gamma), sqrt(ey*bety/gamma), sqrt(ey/bety/gamma), 0.01, 0.002]
###

tmpDir = "_tmp/"
filename = "sigma.dat"

fsigma = open(tmpDir + filename, 'w')

t = twiss2(path + 'doc/FFSexample/assets/ffs.twiss')

print "SLOW: Takes a while to run. Lower the upperbound of the range to 7 for faster results"

for o in range(1,8):

    print "* Order " + str(o)

    # Suppress printing
    old_stdout = sys.stdout
    sys.stdout = DummyFile()
    try:
      m = Map2(t,order=o)
    finally:
      sys.stdout = old_stdout

    mm = Map2(filename=path + 'doc/FFSexample/assets/fort.18',order=o)

    fsigma.write(str(o))

    for i in v:
        s0 = sqrt(mm.sigma(i,sigmaFFS).real)
        sigma = abs(sqrt(m.sigma(i,sigmaFFS).real) - s0)/s0
        fsigma.write("\t" + str(sigma))

    fsigma.write("\n")
