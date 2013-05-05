import sys
from string import split
from operator import *
from numpy import identity, matrix

from definitions import *
import metaclass2
from transport import *

from pytpsa import pol, polmap

from math import *
from multiprocessing import Process
from ctypes import *
cdll.LoadLibrary("../../libs/boost_1_53_0/libboost_python.so.1.53.0")
import mapbeamline_wrapper

################
def gammln(xx):
###############
  g = [0.57236494292474305, 0.0, -0.12078223763524987, -4.4408920985006262e-16, 0.28468287047291829, 0.69314718055994429, 1.2009736023470738, 1.7917594692280547, 2.4537365708424441, 3.1780538303479453, 3.9578139676187165, 4.787491742782044, 5.6625620598571462, 6.5792512120101181, 7.5343642367587762, 8.5251613610654982, 9.5492672573011443, 10.604602902745485, 11.689333420797617, 12.801827480081961, 13.940625219404433, 15.104412573076393, 16.292000476568372, 17.502307845875293, 18.734347511938164, 19.987214495663956, 21.260076156247152, 22.552163853126299, 23.862765841692411, 25.191221182742492, 26.536914491119941, 27.899271383845765, 29.277754515046258, 30.671860106086712, 32.081114895954009, 33.505073450144195, 34.943315776884795, 36.395445208041721, 37.861086508970466, 39.339884187209584]
  return g[int(xx / 0.5 - 1)]


#################
def gammlnGOOD(xx):
#################
  cof = [76.18009172947146, -86.50532032941677, 24.01409824083091, -1.231739572450155, 0.1208650973866179e-2, -0.5395239384953e-5]
  y = x = xx
  tmp = x + 5.5
  tmp -= (x + 0.5) * log(tmp)
  ser = 1.000000000190015
  for c in cof:
    y = y + 1
    ser += c / y
  return - tmp + log(2.5066282746310005 * ser / x)


#########################
class Map2(polmap, dct):
#########################
  '''
  MAP2 coefficients from madx-PTC output

  :param int order: calculate map up to this order
  :param string filename: input filename
  :param twiss2 t: the twiss2 object
  :param twiss2 terr: the twiss object that contains the errors assigned to elements from MADX (e.g. misaligments).

  Either filename is used or twiss object to construct a Map2
  '''

  def __init__(self, *args, **kwargs):
    #if len(args) == 1 and isinstance(args[0], metaclass2.twiss2):
    #  self.fromTwiss(args[0], **kwargs)
    if len(kwargs) == 3 and isinstance(args[0], metaclass2.twiss2):
      self.fromTwissObject(args[0], **kwargs)
    elif len(kwargs) == 4:
      self.fromTwissFile(*args, **kwargs)
    else:
      self.fromFort(*args, **kwargs)

  def fromTwissFile(self, filename, filenameerr=None, order=6, nbProc=1): 
    if filenameerr is None:
      _s = mapbeamline_wrapper.constructMapFromTwissFile(filename, order, nbProc)
    else:
      _s = mapbeamline_wrapper.constructMapFromTwissFileWithErr(filename, filenameerr, order, nbProc) 
    s =  _s.split("|")
    fdct = {}
    for i in range(0, len(s) - 1, 2):
      fdct[str(s[i])] = pol(s[i + 1].strip(" \t\n()[]"), order=order)
    self.update(fdct)
    self.reorder(XYZD)
    
  def fromTwissObject(self, t, terr=None, order=6, nbProc=1):
    if terr is None:
      _s = mapbeamline_wrapper.constructMapFromTwissObject(t, order, nbProc)
    else:
      _s = mapbeamline_wrapper.constructMapFromTwissObjectWithErr(t, terr, order, nbProc) 
    s =  _s.split("|")
    fdct = {}
    for i in range(0, len(s) - 1, 2):
      fdct[str(s[i])] = pol(s[i + 1].strip(" \t\n()[]"), order=order)
    self.update(fdct)
    self.reorder(XYZD)

  ## Twiss
  def fromTwiss(self, t, terr=None, order=6):
    R = generateDefaultMap(order=order)
    U = generateDefaultMatrix(order=order)
    for i in xrange(len(t.elems)):
      e = t.elems[i]
      try:
        mtr = metaclass2.matrixForElement(e, order)
        if mtr == None:
          mp = metaclass2.mapForElement(e, order)
        else:
          M = mtr * U
          mp = metaclass2.matrixToMap(M, XYZD)

        # Apply misalignments here if any
        if terr is not None:
          if isinstance(terr, metaclass2.twiss2):
            dx = terr.elems[i].DX
            dy = terr.elems[i].DY
            if dx != 0:
              mp = mp(x=X+dx)
            if dy != 0:
              mp = mp(y=Y+dy)
          else:
            raise TypeError("The 'align' attribute has to be of the type 'twiss2'.")

        # Combine the map with everything else
        R = mp * R
      except Exception:
        print "No implementation for element: ", e.NAME, e.KEYWORD

    for k in XYZD:
      self[k] = R[k]

    # Reorder the variables so that they are always in the same order
    # This is important for comparision operations but also for all
    # the other methods
    self.reorder(XYZD)

  ## fort.18
  def fromFort(self, order=6, filename='fort.18', nbProc=1):
    ietall = -1

    dct = {}
    fdct = {}

    for line in open(filename):
      sline = split(line)
      l = len(sline)

      if l == 8 or l == 9:
        coef = float(sline[1])
        if l == 8:
          a = [int(sline[3]), int(sline[4]), int(sline[5]), int(sline[6]), int(sline[7])]
        else:
          a = [int(sline[3]), int(sline[4]), int(sline[5]), int(sline[6]), int(sline[7]), int(sline[8])]
        if sum(a) <= order:
          dct[tuple(a)] = coef

      if "etall" in line:
        if ietall >= 0:
          p = pol(order=order)
          p.fromdict(dct, XYZD)
          fdct[XYZD[ietall]] = p
          dct = {}
        ietall += 1

    p = pol()
    p.fromdict(dct, XYZD)
    fdct[XYZD[ietall]] = p
    self.update(fdct)

  def offset(self, xory, sig, gaussianDelta=False):
    '''
    Calculate the beam offset

    :param string xory: Which dimension to calculate for (x,y,px, or py)
    :param list sig: Initial size of beam (sigma) for [x,px,y,py,d,s]
    :param boolean gaussianDelta: Use gaussian energy delta or not

    :return: the offset in the desired dimension
    '''
    sx = 0
    if gaussianDelta:
      if 'x' in xory:
        xory = XYZD[0]
      else:
        xory = XYZD[2]
    for ind, coeff in self[xory].iteritems():
      if all(n % 2 == 0 for n in ind):
        sigmaprod = self.__sigma(ind, sig, gaussianDelta, dv=2)
        if sigmaprod > 0:
          Gammasumln = self.__gamma(ind, gaussianDelta)
          factor = self.__factor(ind, gaussianDelta)
          sx += coeff * factor * exp(Gammasumln) * sigmaprod
    return sx

  def sigma(self, xory, sig, gaussianDelta=False):
    '''
    Calculate the beam size in sigma.

    :param string xory: Which coordinate to calculate for (x,y,px, or py)
    :param list sig: Initial size of beam (sigma) for [x,px,y,py,d,s]
    :param boolean gaussianDelta: Use gaussian energy delta or not

    :return: sigma for xory dimension
    '''
    sx = 0
    for ind1, coeff1 in self[xory].iteritems():
      for ind2, coeff2 in self[xory].iteritems():
        if ind1 >= ind2:
          countfactor = 2.0
          if ind1 == ind2:
            countfactor = 1.0
          ind = [sum(a) for a in zip(ind1, ind2)]
          if all(n % 2 == 0 for n in ind):
            sigmaprod = self.__sigma(ind, sig, gaussianDelta)
            if sigmaprod > 0:
              Gammasumln = self.__gamma(ind, gaussianDelta)
              factor = countfactor * self.__factor(ind, gaussianDelta)
              sx += coeff1 * coeff2 * factor * exp(Gammasumln) * sigmaprod
    return sx

  def comp(self, m, v=None):
    '''
    Compares two maps and returns the chi2. This comparison is only made on the elements for
    which d is 0. It takes two parameters: m for the second map to compare it with and v which
    is a list of the dimensions to be compared, it defaults to the values in XYZD (see definitions.py).

    :param Map2 m: another map
    :param list str v: list of variables used to compare the maps by default it's XYZD

    :return: chi2
    '''
    chi2 = 0
    if v == None:
      v = XYZD[:-2]
    for f in v:
      p1 = self[f]
      p2 = m[f]
      if len(p1.keys()[0]) == len(p2.keys()[0]) and\
         p1.order == p2.order and\
         p1.vars == p2.vars:
        for k, v in p1.iteritems():
          if k[4] == 0:
            chi2 += (v - p2.get(k, 0)) ** 2
      else:
        print "These maps are not comparable."
        return -1
    return chi2

  def compc(self, m, v=None):
    '''
    Compares two maps and returns the chi2. It takes two parameters: m for the second map to
    compare it with and v which is a list of the dimensions to be compared, it defaults to the
    values in XYZD (see definitions.py).


    :param Map2 m: another map
    :param list str v: list of variables used to compare the maps by default it's XYZD

    :return: chi2
    '''
    chi2 = 0
    if v == None:
      v = XYZD[-2]
    for f in v:
      p1 = self[f]
      p2 = m[f]
      if len(p1.keys()[0]) == len(p2.keys()[0]) and\
         p1.order == p2.order and\
         p1.vars == p2.vars:
        for k, v in p1.iteritems():
          chi2 += (v - p2.get(k, 0)) ** 2
      else:
        print "These maps are not comparable."
        return -1
    return chi2

  #Correlation from mapclass.py originally
  def correlation(self, v1, v2, sig, gaussianDelta=False):
    '''
    It calculates the correlation between two dimensions for some given initial sigmas (parameter
    i). Alternatively it can be set to assume a gaussian distribution of the particles.

    :param str v1: name of the first dimension to correlate
    :param str v2: name of the second dimension to correlate
    :param list sig: Initial size of beam (sigma) for [x,px,y,py,d,s]
    :param boolean gaussianDelta: Use gaussian energy delta or not

    :return: correlation
    '''
    sx = 0
    for ind1, coeff1 in self[v1].iteritems():
      for ind2, coeff2 in self[v2].iteritems():
        if ind1 >= ind2 or gaussianDelta:
          countfactor = 2.0
          if ind1 == ind2 or gaussianDelta:
            countfactor = 1.0
          ind = [sum(a) for a in zip(ind1, ind2)]
          if all(n % 2 == 0 for n in ind):
            sigmaprod = self.__sigma(ind, sig, gaussianDelta)
            if sigmaprod > 0:
              Gammasumln = self.__gamma(ind, gaussianDelta)
              factor = countfactor * self.__factor(ind, gaussianDelta)
              sx += coeff1 * coeff2 * factor * exp(Gammasumln) * sigmaprod
    return sx

  #Correlation3 from mapclass.GaussianDelta.py
  def correlation3(self, v1, v2, v3, sig, gaussianDelta=False):
    '''
    It calculates the correlation between two dimensions for some given initial sigmas (parameter
    i). Alternatively it can be set to assume a gaussian distribution of the particles.

    :param str v1: name of the first dimension to correlate
    :param str v2: name of the second dimension to correlate
    :param str v3: name of the third dimension to correlate
    :param list sig: Initial size of beam (sigma) for [x,px,y,py,d,s]
    :param boolean gaussianDelta: Use gaussian energy delta or not

    :return: correlation
    '''
    sx = 0
    for ind1, coeff1 in self[v1].iteritems():
      for ind2, coeff2 in self[v2].iteritems():
        for ind3, coeff3 in self[v3].iteritems():
          countfactor = 1.0
          ind = [sum(a) for a in zip(ind1, ind2, ind3)]
          if all(n % 2 == 0 for n in ind):
            sigmaprod = self.__sigma(ind, sig, gaussianDelta)
            if sigmaprod > 0:
              Gammasumln = self.__gamma(ind, gaussianDelta)
              factor = countfactor * self.__factor(ind, gaussianDelta)
              sx += coeff1 * coeff2 * coeff3 * factor * exp(Gammasumln) * sigmaprod
    return sx

  def generatelist(self, xory, sig, gaussianDelta=False):
    '''
    Provides a list of the largest contributions to the sigmas assuming a Gaussian distribution
    or a uniform one.

    :param str xory: name of the dimension
    :param list sig: Initial size of beam (sigma) for [x,px,y,py,d,s]
    :param boolean gaussianDelta: Use gaussian energy delta or not

    :return: list of contributions
    '''
    sx = 0
    l = []
    if gaussianDelta:
      if 'x' in xory:
        xory = XYZD[0]
      else:
        xory = XYZD[2]
    for ind1, coeff1 in self[xory].iteritems():
      for ind2, coeff2 in self[xory].iteritems():
        if ind1 >= ind2:
          countfactor = 2.0
          if ind1 == ind2:
            countfactor = 1.0
          ind = [sum(a) for a in zip(ind1, ind2)]
          if all(n % 2 == 0 for n in ind):
            sigmaprod = self.__sigma(ind, sig, gaussianDelta)
            if sigmaprod > 0:
              Gammasumln = self.__gamma(ind, gaussianDelta)
              factor = countfactor * self.__factor(ind, gaussianDelta)
              sxt = coeff1 * coeff2 * factor * exp(Gammasumln) * sigmaprod
              l.append([-abs(sxt), sxt] + list(ind1) + list(ind2))
    l.sort()
    return l

  #Auxiliary functions (private)
  def __sigma(self, ind, i, gaussianDelta, dv=1):
    if gaussianDelta:
#      sigmaprod = pow(i[0], ind[0])*pow(i[1], ind[1])*pow(i[2], ind[2])*pow(i[3], ind[3])*pow(i[4], ind[4])
      sigmaprod = reduce(mul, map(pow, i, ind))
    else:
#      sigmaprod = pow(i[0], ind[0])*pow(i[1], ind[1])*pow(i[2], ind[2])*pow(i[3], ind[3])*pow(i[4]/2., ind[4])
      qq = 1
      if len(ind) > 5:
        qq = pow(i[5] / dv, ind[5])
      sigmaprod = reduce(mul, map(pow, i[:4], ind[:4])) * pow(i[4] / 2., ind[4]) * qq
    return sigmaprod

  def __gamma(self, ind, gaussianDelta):
    if gaussianDelta:
#      Gammasumln = gammln(0.5+ind[0]/2.)+gammln(0.5+ind[1]/2.)+gammln(0.5+ind[2]/2.)+gammln(0.5+ind[3]/2.)+gammln(0.5+ind[4]/2.)
      Gammasumln = reduce(add, map(lambda x: gammln(0.5 + x / 2.), ind))
    else:
      indl = list(ind)
      del indl[4]
      Gammasumln = reduce(add, map(lambda x: gammln(0.5 + x / 2.), indl))
    return Gammasumln

  def __factor(self, ind, gaussianDelta, poten=2.):
    l = len(ind)
    if gaussianDelta:
        if l == 5:
          poten = 2.5
        if l == 6:
          poten = 3.0
        factor = pow(2, sum(ind) / 2.) / pow(pi, poten)
    else:
      if l > 5:
        poten = 2.5
      factor = pow(2, sum(ind[:4] + ind[5:]) / 2) / pow(pi, poten) / (ind[4] + 1.0)
    return factor
