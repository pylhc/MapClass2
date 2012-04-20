from string import split
from math import *

from pytpsa import pol, polmap


################
def gammln(xx):
###############
  g=[0.57236494292474305,0.0,-0.12078223763524987,-4.4408920985006262e-16,0.28468287047291829,0.69314718055994429,1.2009736023470738,1.7917594692280547,2.4537365708424441,3.1780538303479453,3.9578139676187165,4.787491742782044,5.6625620598571462,6.5792512120101181,7.5343642367587762,8.5251613610654982,9.5492672573011443,10.604602902745485,11.689333420797617,12.801827480081961,13.940625219404433,15.104412573076393,16.292000476568372,17.502307845875293,18.734347511938164,19.987214495663956,21.260076156247152,22.552163853126299,23.862765841692411,25.191221182742492,26.536914491119941,27.899271383845765,29.277754515046258,30.671860106086712,32.081114895954009,33.505073450144195,34.943315776884795,36.395445208041721,37.861086508970466,39.339884187209584]
  return g[int(xx/0.5-1)]



#################
def gammlnGOOD( xx):
#################
  cof=[76.18009172947146,-86.50532032941677,24.01409824083091,-1.231739572450155,0.1208650973866179e-2,-0.5395239384953e-5]
  y=x=xx
  tmp=x+5.5
  tmp -= (x+0.5)*log(tmp)
  ser=1.000000000190015;
  for c in cof:
    y=y+1
    ser += c/y
  return -tmp+log(2.5066282746310005*ser/x)



#########################
class Map2(polmap):
#########################
  '''
  MAP coefficients from madx-PTC output

  :param int order: Calculate map up to this order
  :param string filename: Input filename
  '''

  def __init__(self, order=6, filename='fort.18'):
    ietall=-1
    
    xyzd=['x', 'px', 'y', 'py', 'd','s']
    fxyzd=['fx', 'fpx', 'fy', 'fpy', 'fd', 'fs']
    dct={}
    fdct={}

    for line in open(filename):
      sline=split(line)
      l=len(sline)
      # Filling sline with 0's at the end, that way it doesn't matter if we read a file with 5 vars 
      # or 6 vars. If it's 5 vars the order for 's' will always be 0 and give the right results
      sline+=[0,0]

      if l == 8 or l == 9:
        coef=float(sline[1])
        a=[int(sline[3]), int(sline[4]), int(sline[5]), int(sline[6]), int(sline[7]), int(sline[8])]
        if sum(a) <= order:
          dct[tuple(a)]=coef

      if "etall" in line:
        if ietall >= 0:
          p=pol()
          p.fromdict(dct,xyzd)
          fdct[fxyzd[ietall]]=p
          dct={}
        ietall+=1

    p=pol()
    p.fromdict(dct,xyzd)
    fdct[fxyzd[ietall]]=p
    self.update(fdct)


  def offset(self, xory, i, gaussianDelta=False):
    '''
    Calculate the beam offset

    :param string xory: Which coordinate to calculate for (x,y,px, or py)
    :param list i: Size of beam in sigma [x,px,y,py]
    '''
    sx=0
    for ind,coeff in self[xory].iteritems():
      if all(n % 2 == 0 for n in ind):
        sigmaprod=self.__sigma(ind, i, gaussianDelta)
        if sigmaprod > 0:
          Gammasumln=self.__gamma(ind, gaussianDelta)
          factor=self.__factor(ind, gaussianDelta)
          sx+=coeff*factor*exp(Gammasumln)*sigmaprod
    return sx


  def sigma(self, xory, i, gaussianDelta=False):
    '''
    Calculate the beam size in sigma.

    :param string xory: Which coordinate to calculate for (x,y,px, or py)
    :param list i: Size of beam in sigma [x,px,y,py]
    '''
    sx=0
    for ind1,coeff1 in self[xory].iteritems():
      for ind2,coeff2 in self[xory].iteritems():
        if ind1 >= ind2:
          countfactor=2.0
          if ind1 == ind2:
            countfactor=1.0
          ind=[sum(a) for a in zip(ind1, ind2)]
          if all(n % 2 == 0 for n in ind):
            sigmaprod=self.__sigma(ind, i, gaussianDelta)
            if sigmaprod > 0:
              Gammasumln=self.__gamma(ind, gaussianDelta)
              factor=countfactor*self.__factor(ind, gaussianDelta)
              sx+=coeff1*coeff2*factor*exp(Gammasumln)*sigmaprod
    return sx

  #Auxiliary functions (private)
  def __sigma(self, ind, i, gaussianDelta):
    if (gaussianDelta):
      sigmaprod = pow(i[0], ind[0])*pow(i[1], ind[1])*pow(i[2], ind[2])*pow(i[3], ind[3])*pow(i[4], ind[4])
    else:
      sigmaprod = pow(i[0], ind[0])*pow(i[1], ind[1])*pow(i[2], ind[2])*pow(i[3], ind[3])*pow(i[4]/2., ind[4])
    return sigmaprod


  def __gamma(self, ind, gaussianDelta):
    if (gaussianDelta):
      Gammasumln = gammln(0.5+ind[0]/2.)+gammln(0.5+ind[1]/2.)+gammln(0.5+ind[2]/2.)+gammln(0.5+ind[3]/2.)+gammln(0.5+ind[4]/2.)
    else:
      Gammasumln = gammln(0.5+ind[0]/2.)+gammln(0.5+ind[1]/2.)+gammln(0.5+ind[2]/2.)+gammln(0.5+ind[3]/2.)
    return Gammasumln

  def __factor(self, ind, gaussianDelta):
    if (gaussianDelta):
      factor = pow(2, sum(ind)/2.)/pow(pi, 2.5)
    else:
      #FIXME: This works for now but it's not right!!!
      factor = pow(2, sum(ind[:-2])/2.)/pow(pi, 2.)/(ind[-2]+1)
    return factor
