import sys

sys.path.append("libs")

import os
from math import *

try:
    from collections import OrderedDict
except ImportError:
    print "Python 2.7+ OrderedDict collection not available"
    try:
        from ordereddict import OrderedDict
        print "Using backported OrderedDict implementation"
    except ImportError:
        print "Backported OrderedDict implementation not available"

import itertools
import argparse

import unittest
if sys.version_info[0] == 2 and sys.version_info[1] < 7:
    import unittest2
else:
    unittest2 = unittest
    # probably won't work w/ Python 3.0 / 3.1, tested on Python 3.2 and 2.7
    import unittest

from definitions import *
from mapclass import Map2

sys.path.append("old")


# Make sets of strings to compare them because oldMapClass might
# return less elements than MapClass2
#
# We use strings and not floats because the accuracy for floats is not
# very good and there are some discrepancies between the two versions
# after the 10th decimal

class TestExtended:

  xyzd = XYZD[:-2]

  sf = "%.8e"

  compare = False
  correlation = False
  gaussian = False

  ### Example values
  betx = 66.14532014
  bety = 17.92472388
  gamma = 3e6
  ex = 68e-8
  ey = 2e-8
  sigmaFFS = [sqrt(ex*betx/gamma), sqrt(ex/betx/gamma), sqrt(ey*bety/gamma), sqrt(ey/bety/gamma), 0.01]
  ###

  def assertAlmostEq(self, a, b):
    self.assertEqual(self.sf % a, self.sf % b)

  def testEqual(self):
    res1 = set([self.sf % e[()] for e in self.m(self.vals).values()])
    res2 = set([self.sf % e for e in self.mm.f(self.vals.values())])
    self.assertTrue(res1.issuperset(res2))

  def testSigma(self):
    for i in self.xyzd:
      self.assertAlmostEq(self.m.sigma(i,self.sigmaFFS,self.gaussian), self.mm.sigma(i,self.sigmaFFS))

  def testOffset(self):
    for i in self.xyzd:
      self.assertAlmostEq(self.m.offset(i,self.sigmaFFS,self.gaussian), self.mm.offset(i,self.sigmaFFS))

  def testCorrelation(self):
    if self.correlation:
      for [v1, v2] in itertools.combinations(self.xyzd, 2):
        self.assertAlmostEq(self.m.correlation(v1,v2,self.sigmaFFS,self.gaussian), self.mm.correlation(v1,v2,self.sigmaFFS))

  @unittest2.skipUnless('-s' in sys.argv, "Slow test")
  def testCorrelation3(self):
    if self.correlation:
      for [v1, v2, v3] in itertools.combinations(self.xyzd, 3):
        self.assertAlmostEq(self.m.correlation3(v1,v2,v3,self.sigmaFFS,self.gaussian), self.mm.correlation3(v1,v2,v3,self.sigmaFFS))

  def testComp(self):
    if self.compare:
      self.assertAlmostEq(self.m.comp(self.m2,v=self.xyzd), self.mm.comp(self.mm2))

  def testCompc(self):
    if self.compare:
      self.assertAlmostEq(self.m.compc(self.m2,v=self.xyzd), self.mm.compc(self.mm2))

  def testGenList(self):
    for v in self.xyzd:
      self.mm.generatelist(v, self.sigmaFFS)
      for (l1,l2) in zip(self.m.generatelist(v, self.sigmaFFS, self.gaussian), getattr(self.mm, 'list' + v)):
        self.assertAlmostEq(l1[0],l2[0])
        # Because the order in which the indices are stored isn't the
        # same the result isn't equal but equivalent.
        # E.g.
        # l1[2:] = (1,2,3,0,1,0) equal l2[2:] = (1,2,3,0,1,0)
        # l1[2:] = (1,2,3,0,1,0) equivalent l2[2:] = (0,1,0,1,2,3)
        ind1 = l1[2:]
        ind2 = l2[2:]
        l = len(ind1)/2
        self.assertTrue((ind1[:l] == ind2[:l] or ind1[:l] == ind2[l:]) and (ind1[l:] == ind2[:l] or ind1[l:] == ind2[l:]))


################
## Test cases ##
################

class Test5var6order(unittest2.TestCase, TestExtended):

  def setUp(self):
    from mapclass25 import Map

    o = 6
    f = 'assets/fort.18'
    ff = 'assets/5fort.18'
    self.compare = True
    self.vals = OrderedDict([('x',self.sigmaFFS[0]),('px',self.sigmaFFS[1]),('y',self.sigmaFFS[2]),('py',self.sigmaFFS[3]),('d',self.sigmaFFS[4]),('s',0)])

    self.m = Map2(order=o, filename=f)
    self.mm = Map(order=o, filename=f)
    self.m2 = Map2(order=o, filename=ff)
    self.mm2 = Map(order=o, filename=ff)


class Test5var6orderGaussian(unittest2.TestCase, TestExtended):

  def setUp(self):
    from mapclassGaussianDelta25 import Map

    o = 6
    f = 'assets/fort.18'
    self.gaussian = True
    self.correlation = True
    self.vals = OrderedDict([('x',self.sigmaFFS[0]),('px',self.sigmaFFS[1]),('y',self.sigmaFFS[2]),('py',self.sigmaFFS[3]),('d',self.sigmaFFS[4]),('s',0)])

    self.m = Map2(order=o, filename=f)
    self.mm = Map(order=o, filename=f)


class Test5var10order(unittest2.TestCase, TestExtended):

  def setUp(self):
    from mapclass25 import Map

    o = 10
    f = 'assets/fort.18'
    ff = 'assets/5fort.18'
    self.compare = True
    self.vals = OrderedDict([('x',self.sigmaFFS[0]),('px',self.sigmaFFS[1]),('y',self.sigmaFFS[2]),('py',self.sigmaFFS[3]),('d',self.sigmaFFS[4]),('s',0)])

    self.m = Map2(order=o, filename=f)
    self.mm = Map(order=o, filename=f)
    self.m2 = Map2(order=o, filename=ff)
    self.mm2 = Map(order=o, filename=ff)


class Test5var10orderGaussian(unittest2.TestCase, TestExtended):

  def setUp(self):
    from mapclassGaussianDelta25 import Map

    o = 10
    f = 'assets/fort.18'
    self.gaussian = True
    self.correlation = True
    self.vals = OrderedDict([('x',self.sigmaFFS[0]),('px',self.sigmaFFS[1]),('y',self.sigmaFFS[2]),('py',self.sigmaFFS[3]),('d',self.sigmaFFS[4]),('s',0)])

    self.m = Map2(order=o, filename=f)
    self.mm = Map(order=o, filename=f)


class Test6var6order(unittest2.TestCase, TestExtended):

  def setUp(self):
    from mapclass25_6var import Map

    self.sigmaFFS = [0.000227, 9.306e-5, 5.02e-4, 4.21e-5, 0.00666, 0.002]

    o = 6
    f = 'assets/6Dfort.18'
    self.vals = OrderedDict([('x',self.sigmaFFS[0]),('px',self.sigmaFFS[1]),('y',self.sigmaFFS[2]),('py',self.sigmaFFS[3]),('d',self.sigmaFFS[4]),('s',self.sigmaFFS[5])])

    self.m = Map2(order=o, filename=f)
    self.mm = Map(order=o, filename=f)

  @unittest2.skip("Unnecessary")
  def testGenList(self):
    # generatelist isn't clear in mapclass25_6var so override this
    # test to avoid running it
    return


class Test6var10order(unittest2.TestCase, TestExtended):

  def setUp(self):
    from mapclass25_6var import Map

    self.sigmaFFS = [0.000227, 9.306e-5, 5.02e-4, 4.21e-5, 0.00666, 0.002]

    o = 10
    f = 'assets/6Dfort.18'
    self.vals = OrderedDict([('x',self.sigmaFFS[0]),('px',self.sigmaFFS[1]),('y',self.sigmaFFS[2]),('py',self.sigmaFFS[3]),('d',self.sigmaFFS[4]),('s',self.sigmaFFS[5])])

    self.m = Map2(order=o, filename=f)
    self.mm = Map(order=o, filename=f)

  @unittest2.skip("Unnecessary")
  def testGenList(self):
    # generatelist isn't clear in mapclass25_6var so override this
    # test to avoid running it
    return


########################
## Twiss
########################

class TwissExtended:

  xyzd = XYZD[:-2]

  def compareFortTwiss(self):
    self.assertTrue(self.m.compc(self.mm, self.xyzd).real < 0.1)

  def testBeta(self):
    for i in range(len(self.t.elems)):
      e = self.t.elems[i]
      if e.L != 0:
          self.assertEqual(e.BETX, self.t.getBeta(i, e.L).BETX)
          self.assertEqual(e.BETY, self.t.getBeta(i, e.L).BETY)
          self.assertEqual(e.ALFX, self.t.getBeta(i, e.L).ALFX)
          self.assertEqual(e.ALFY, self.t.getBeta(i, e.L).ALFY)

  def testDisp(self):
    for i in range(len(self.t.elems)):
      e = self.t.elems[i]
      self.assertEqual(e.DX, self.t.getDisp(i, e.L).DX)
      self.assertEqual(e.DY, self.t.getDisp(i, e.L).DY)
      self.assertEqual(e.DPX, self.t.getDisp(i, e.L).DPX)
      self.assertEqual(e.DPY, self.t.getDisp(i, e.L).DPY)


class TestElems(unittest2.TestCase, TwissExtended):

  def setUp(self):
    from metaclass2 import twiss2

    for root, subFolders, files in os.walk('assets/twiss/'):
      if 'twiss' in files and 'fort.18' in files:
        print "\nTesting: ", root
        twissfile = os.path.join(root, "twiss")
        fortfile = os.path.join(root, "fort.18")
        self.t = twiss2(twissfile)
        self.m = Map2(self.t)
        self.mm = Map2(filename=fortfile)


def twissSuite():
  suite = unittest2.TestSuite()
  suite.addTest(TestElems("compareFortTwiss"))
  suite.addTest(TestElems("testBeta"))
  suite.addTest(TestElems("testDisp"))
  return suite

####################
## Main
####################

if __name__ == '__main__':
  parser = argparse.ArgumentParser(description='Tester for MapClass')
  parser.add_argument('-s', help='Run the slow tests as well', action='store_true', dest='slow')
  parser.add_argument('-t', help="Test the import from Twiss and compare with fort.18. It doesn't run any of the other tests.", action='store_true', dest='twiss')
  args = parser.parse_args()

  if args.twiss:
    runner = unittest2.TextTestRunner(verbosity=2)
    test_suite = twissSuite()
    runner.run(test_suite)
  else:
    unittest2.main(verbosity=2, argv=sys.argv[:1])
