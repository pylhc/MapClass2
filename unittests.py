#!/usr/bin/env python

import sys
from math import *
from collections import OrderedDict
import itertools
import argparse
import unittest

from mapclass import Map2

sys.path.append("old")

# Make sets of strings to compare them because oldMapClass might
# return less elements than MapClass2
#
# We use strings and not floats because the accuracy for floats is not
# very good and there are some discrepancies between the two versions
# after the 10th decimal

class TestExtended:

  fxyzd = ['fx', 'fpx', 'fy', 'fpy']

  sf = "%.10e"

  compare = False
  correlation = False
  gaussian = False

  ###
  betx=66.14532014
  bety=17.92472388
  gamma=3e6
  ex=68e-8
  ey=2e-8
  sigmaFFS=[sqrt(ex*betx/gamma), sqrt(ex/betx/gamma), sqrt(ey*bety/gamma), sqrt(ey/bety/gamma), 0.01]
  ###

  def assertAlmostEq(self,a,b):
    self.assertEqual(self.sf % a, self.sf % b)

  def testEqual(self):
    res1 = set([self.sf % e[()] for e in self.m(self.vals).values()])
    res2 = set([self.sf % e for e in self.mm.f(self.vals.values())])
    self.assertTrue(res1.issuperset(res2))

  def testSigma(self):
    for i in self.fxyzd:
      self.assertAlmostEq(self.m.sigma(i,self.sigmaFFS,self.gaussian), self.mm.sigma(i[1:],self.sigmaFFS))

  def testOffset(self):
    for i in self.fxyzd:
      self.assertAlmostEq(self.m.offset(i,self.sigmaFFS,self.gaussian), self.mm.offset(i[1:],self.sigmaFFS))

  def testCorrelation(self):
    if self.correlation:
      for [v1, v2] in itertools.combinations(self.fxyzd, 2):
        self.assertAlmostEq(self.m.correlation(v1,v2,self.sigmaFFS,self.gaussian), self.mm.correlation(v1[1:],v2[1:],self.sigmaFFS))

  @unittest.skipUnless('-s' in sys.argv, "Slow test")
  def testCorrelation3(self):
    if self.correlation:
      for [v1, v2, v3] in itertools.combinations(self.fxyzd, 3):
        self.assertAlmostEq(self.m.correlation3(v1,v2,v3,self.sigmaFFS,self.gaussian), self.mm.correlation3(v1[1:],v2[1:],v3[1:],self.sigmaFFS))

  def testComp(self):
    if self.compare:
      self.assertAlmostEq(self.m.comp(self.m2,v=self.fxyzd), self.mm.comp(self.mm2))

  def testCompc(self):
    if self.compare:
      self.assertAlmostEq(self.m.compc(self.m2,v=self.fxyzd), self.mm.compc(self.mm2))

################
## Test cases ##
################

class Test5var6order(unittest.TestCase, TestExtended):

  def setUp(self):
    from mapclass25 import Map

    o = 6
    f = 'assets/fort.18'
    ff = 'assets/5fort.18'
    self.compare = True
    self.vals = OrderedDict([('x',1E-19),('px',3E-23),('y',5E-32),('py',6E-12),('d',-1.5E23),('s',0)])

    self.m = Map2(order=o,filename=f)
    self.mm = Map(order=o,filename=f)
    self.m2 = Map2(order=o,filename=ff)
    self.mm2 = Map(order=o,filename=ff)

class Test5var6orderGaussian(unittest.TestCase, TestExtended):

  def setUp(self):
    from mapclassGaussianDelta25 import Map

    o = 6
    f = 'assets/fort.18'
    self.gaussian = True
    self.correlation = True
    self.vals = OrderedDict([('x',1E-19),('px',3E-23),('y',5E-32),('py',6E-12),('d',-1.5E23),('s',0)])

    self.m = Map2(order=o,filename=f)
    self.mm = Map(order=o,filename=f)


if __name__ == '__main__':
  parser = argparse.ArgumentParser(description='Tester for MapClass')
  parser.add_argument('-s', help='Run the slow tests as well',
                    action='store_true', dest='slow')
  parser.parse_args()

  unittest.main(verbosity=2,argv=sys.argv[:1])
