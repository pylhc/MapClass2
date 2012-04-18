from string import split

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
class Map2:
#########################
  '''
  MAP coefficients from madx-PTC output

  :param int order: Calculate map up to this order
  :param string filename: Input filename
  :param boolean gaussianDelta: Use gaussianDelta or not
  '''

  def __init__(self, order=6, filename='fort.18', gaussianDelta=False):
    ietall=0
    self.order=order+1
    self.gaussianDelta=gaussianDelta
    
    xyzd=['x', 'px', 'y', 'py', 'd','s']
    fxyzd=['fx', 'fpx', 'fy', 'fpy', 'fd', 'fs']
    dct={}
    fdct={}

    for line in open(filename):
      sline=split(line)

      if len(sline) == 8:
        coef=float(sline[1])
        a=[int(sline[3]), int(sline[4]), int(sline[5]), int(sline[6]), int(sline[7])]
        if (a[0]+a[1]+a[2]+a[3]+a[4]) <= self.order:
          dct[tuple(a)]=coef

      if "etall" in line:
        p=pol()
        p.fromdict(dct,xyzd)
        fdct[fxyzd[ietall]]=p
        dct={}
        ietall+=1

    self.polmap=polmap(fdct)
