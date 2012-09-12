import sys

sys.path.append('../../')

from metaclass2 import twiss2
from mapclass import Map2
from math import sqrt

######## EXAMPLE SCRIPT ########
# Varies position of all sextupoles in the beamline in 0.1m steps and
# recalculates sigmas each time. If sigma has been reduced, prints the
# position change of the sextupole. For each sextupole the sigmas are
# compared with those of the original twiss object.

# Constants
order = 2

betx = 64.99988501
bety = 17.99971417
gamma = 3e6
ex = 68e-8
ey = 2e-8
sigmaFFS = [sqrt(ex*betx/gamma), sqrt(ex/betx/gamma), sqrt(ey*bety/gamma), sqrt(ey/bety/gamma), 0.01, 0.002]

# Initialising
t0 = twiss2('assets/ffs.twiss')
print "* Stripping the line."
t0 = t0.stripLine()
print "* Merging the elements."
t0 = t0.mergeElems()
print "* Calculating the map."
m0 = Map2(t0, order=order)
print "* Calculating original sigma x ",
sOrig = sqrt(m0.sigma('x', sigmaFFS).real)
print sOrig

# Loop through the elements in the stripped and merged beamline
print "* Starting main loop"
for i in range(1, len(t0.elems)-1):
    # Define current, previous and subsequent elements and continue if
    # i is a sextupole surrounded by drifts
    e = t0.elems[i]
    prev = t0.elems[i-1]
    nxt = t0.elems[i+1]
    if e.KEYWORD == "SEXTUPOLE" and prev.KEYWORD == "DRIFT" and nxt.KEYWORD == "DRIFT":
        print ""
        print "!!! Sextupole found with element No. %d. Looking for optimisations" % i,
        # Define lower and upper bounds as lengths of surrounding drifts
        lower = prev.L * -1
        upper = nxt.L
        # For each sextupole the original sigma is used for comparison
        # dPosIdeal is the change in position resulting in the
        # smallest sigma found for the current sextupole
        s0 = sOrig
        dPosIdeal = 0
        # Vary position of sextupole from lower bound to upper bound
        # in increments of 0.1m. Recalculates twiss, map and sigmas
        # having altered the element's position
        while lower < upper:
            sys.stdout.write(".")
            sys.stdout.flush()
            t = t0.alterElem(i, dPos=lower)
            if t is None: break # Stops if alterElem finds an issue
            m = Map2(t,order=order)
            s = sqrt(m.sigma('x', sigmaFFS).real)
            # If the new sigma is smaller than the previous one, make
            # them equal to compare with the next position alteration
            # and store the position change.
            if s < s0:
                s0 = s
                dPosIdeal = lower
            lower += 0.1

        print ""
        if dPosIdeal != 0:
            print "Original sigma x: ", sOrig
            print "Optimised sigma x: ", s0
            print "Sigma x reduced by: ", sOrig - s0
            print "This reduction occurred at a position change of: ", dPosIdeal
        else:
            print "This sextupole could not be optimised"
