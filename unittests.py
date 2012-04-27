import sys
from math import *
from collections import OrderedDict
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

  gaussian = False

  def testEqual(self):
    res1 = set([self.sf % e[()] for e in self.m(self.vals).values()])
    res2 = set([self.sf % e for e in self.mm.f(self.vals.values())])
    self.assertTrue(res1.issuperset(res2))

  def testSigma(self):
    gaussian=self.gaussian
    betx=66.14532014
    bety=17.92472388
    gamma=3e6
    ex=68e-8
    ey=2e-8
    sigmaFFS=[sqrt(ex*betx/gamma), sqrt(ex/betx/gamma), sqrt(ey*bety/gamma), sqrt(ey/bety/gamma), 0.01]
    for i in self.fxyzd:
      self.assertEqual(self.sf % self.m.sigma(i,sigmaFFS,gaussian),self.sf % self.mm.sigma(i[1:],sigmaFFS))

  def testOffset(self):
    gaussian=self.gaussian
    betx=66.14532014
    bety=17.92472388
    gamma=3e6
    ex=68e-8
    ey=2e-8
    sigmaFFS=[sqrt(ex*betx/gamma), sqrt(ex/betx/gamma), sqrt(ey*bety/gamma), sqrt(ey/bety/gamma), 0.01]
    for i in self.fxyzd:
      print i
      self.assertEqual(self.sf % self.m.offset(i,sigmaFFS,gaussian),self.sf % self.mm.offset(i[1:],sigmaFFS))



class Test5var6order(unittest.TestCase, TestExtended):
  
  def setUp(self):
    from mapclass25 import Map

    self.o = 6
    self.f = 'assets/fort.18'
    self.vals = OrderedDict([('x',1E-19),('px',3E-23),('y',5E-32),('py',6E-12),('d',-1.5E23),('s',0)])

    self.m = Map2(order=self.o,filename=self.f)
    self.mm = Map(order=self.o,filename=self.f)

class Test5var6orderGaussian(unittest.TestCase, TestExtended):

  def setUp(self):
    from mapclassGaussianDelta25 import Map

    self.o = 6
    self.f = 'assets/fort.18'
    self.gaussian = True
    self.vals = OrderedDict([('x',1E-19),('px',3E-23),('y',5E-32),('py',6E-12),('d',-1.5E23),('s',0)])

    self.m = Map2(order=self.o,filename=self.f)
    self.mm = Map(order=self.o,filename=self.f)


if __name__ == '__main__':
  unittest.main()
