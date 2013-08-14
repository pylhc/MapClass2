import sys
from os import system
from ctypes import *
cdll.LoadLibrary("../../libs/boost_1_53_0/libboost_python.so.1.53.0")
sys.path.append('../../')

from metaclass2 import twiss2
from mapclass import Map2
from math import sqrt
from matplotlib import pyplot as plt
from numpy import arange


######## EXAMPLE SCRIPT ########
# asks for user input to search for optimised strength given sextupole, range, and increment

# Constants
order = 3
gaussian = False
betx = 64.99988501
bety = 17.99971417
gamma = 3e6
ex = 68e-8
ey = 2e-8


# Initialising
t0 = twiss2('assets/ffs.twiss')
print "* Stripping the line."
t0 = t0.stripLine()
print "* Merging the elements."
t0 = t0.mergeElems()
print "* Calculating the map."
# Determine number of sextupole to be analyzed
T = 1
while T == 1:
    var = raw_input("* Enter the element number of sextupole to be analyzed: ")
    try:
        i = int(var)
        if i >= len(t0.elems):
            print "* Element numbers are less than", len(t0.elems),"."
        elif i < 0:
            print "* Element numbers are positive."
        else:
            print "* You have chosen element", i
            if t0.elems[i].KEYWORD != "SEXTUPOLE":
                print "* This element is not a sextupole."
            else:
                sxtnum = i
                T = 0    
    except ValueError:
        print "* This is not a valid element number. Input must be an integer."
while T == 0:
    beg = float(raw_input("* Where would you like to begin the search range? "))
    end = float(raw_input("* Where would you like to end the search range? "))
    step = float(raw_input("* What step size would you like to use? "))
    try:
        posrange = arange(beg, end, step)
        T = 1
    except ValueError:
        print "* Range or step size are invalid."
        T = 0
        
#Loop through values of dPos for given element

chilist = [0]*len(posrange)
print "* Calculating chi-squared for modified position of element", sxtnum,"."
for i in range(len(posrange)):
    sys.stdout.write(".")
    sys.stdout.flush()
    t = t0.alterElem(sxtnum, dPos = posrange[i])
    if t is None:
        break
    sigmaFFS = [sqrt(ex*betx/gamma), sqrt(ex/betx/gamma), sqrt(ey*bety/gamma), sqrt(ey/bety/gamma), 0.01, 0.002]
    m0 = Map2(t, terr=None, order = order, nbProc=2)

    h0 = 40e-9 #sigmaFFS[0]
    j0 = 1e-9  #sigmaFFS[2]
    h = sqrt(m0.sigma('x', sigmaFFS, gaussian).real)
    j = sqrt(m0.sigma('y', sigmaFFS, gaussian).real)
    chi = (h-h0)**2/40**2 + (j-j0)**2
    #chi = (h/40)**2 + (j)**2
    
    chilist[i] = chi
chimin = min(chilist)
plt.plot(posrange, chilist)
plt.ylabel('chi-squared')
plt.xlabel('position')
plt.title('Position optimisation for element %d' % sxtnum)
for i in range(len(posrange)):
    if chilist[i] == chimin:
        print "* chi-square minimized to", chimin, "at pos = ", posrange[i]
plt.savefig('chiposrange_elem_%d.png' % sxtnum)
plt.show()
