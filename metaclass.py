from collections import namedtuple
from numpy import identity
from copy import copy

from transport import *
import mapclass
from pytpsa import polmap
import math
from integration import simpson

#########################
# Auxiliar elements
#########################
class dct(dict):
  def __getattr__(self, attr):
    try:
      return self[attr]
    except:
      raise AttributeError("%r object has no attribute %r" % (type(self).__name__, attr))

#########################
class twiss(dict):
#########################
  "Twiss parameters from madx output (with free choice of select items)"

  def __init__(self, filename='twiss'):
    self.elems = []
    # Markers
    self.markers = []

    f = open(filename, 'r')

    for line in f:
      # if ("@ " in line and "%le" in line) :    # FIX to take DPP %s
      if "@ " not in line and "@" in line:
        line = replace(line, "@" , "@ ")

      splt = line.split()
      if "@ " in line and "%" in line and "s" not in splt[2]:
        label = splt[1]
        try:
            self[label] = float(splt[3])
        except:
            print "Problem parsing:", line
            print "Going to be parsed as string"
            try:
              self[label] = splt[3].strip('"')
            except:
              print "Problem persits, let's ignore it!"

      elif "@ " in line and "s" in splt[2]:
        label = splt[1].strip(":")
        self[label] = splt[3].strip('"')

      if "* " in line or "*\t" in line:
        labels = splt[1:]
        print "labels ",len(labels)

      if "$ " in line or "$\t" in line:
        types = splt[1:]

      if "@" not in line and "*" not in line and "%" not in line:
        vals = []
        for j in range(0,len(labels)):
          if "%hd" in types[j]:
            vals.append(int(splt[j]))
          if "%le" in types[j]:
            vals.append(float(splt[j]))
          if "s" in types[j]:
            vals.append(splt[j].strip('"'))
        e = dct(zip(labels,vals))
        if "$" in line:
          self.markers.append(e)
        else:
          self.elems.append(e)

    f.close()
    try:
      labels
      types
    except:
      print "From Metaclass: Bad format or empty file ", filename
      print "Leaving Metaclass"
      exit()

  ## BETA FUNCTION: Find beta, alpha and gamma functions at location s in element nE
  # nE = number of element of interest, s = location of interest within element nE
  def getBeta(self,nE,s):
    # Get initial beta, alpha and gamma at end of previous element
    # If nE is first element, use start marker as initial; if nE any other element, use previous element
    if nE != 0:
      prevE = self.elems[nE-1]
    else:
      prevE = self.markers[0]

    betX0 = prevE.BETX
    alfX0 = prevE.ALFX
    gamX0 = (1+alfX0**2)/betX0

    betY0 = prevE.BETY
    alfY0 = prevE.ALFY
    gamY0 = (1+alfY0**2)/betY0

    if s == 0:
      return  dct([ ('BETX',betX0),
                    ('ALFX',alfX0),
                    ('GAMX',gamX0),
                    ('BETY',betY0),
                    ('ALFY',alfY0),
                    ('GAMY',gamY0) ])

    paraX0 = mtrx([ [betX0],
                    [alfX0],
                    [gamX0] ])

    paraY0 = mtrx([ [betY0],
                    [alfY0],
                    [gamY0] ])

    paraInitial = [paraX0, paraY0]

    # Copy element of interest and change its length to location of interest, s
    # Calculate transport matrix for this element assuming D=0
    # If nE is not DRIFT, QUADRUPOLE or DIPOLE, change element to DRIFT and recalculate transport matrix
    e = dct(self.elems[nE])
    e['L'] = s
    eTransport = matrixForElement(e, 6)    # NOTE: Take order 6
    if eTransport == None:
      e['KEYWORD']="DRIFT"
      eTransport = matrixForElement(e, 6)
    eTransport = eTransport(d=0)

    # Extract required elements from transport matrix to form transform matrix
    # i=0 for horizontal x-direction; i=1 for vertical y-direction
    # a=C, b=S, c=C', d=S' according to standard notation
    para = []
    for i in [0, 1]:
      j = 2*i
      a = eTransport.item((j,j))    # Equivalent to item((0,0)) for x-direction; item((2,2)) for y-direction
      b = eTransport.item((j,j+1))
      c = eTransport.item((j+1,j))
      d = eTransport.item((j+1,j+1))
      twissTransform = mtrx([ [a**2, -2*a*b, b**2],  # Create transform matrix
                              [-a*c, b*c + a*d, -b*d],
                              [c**2, -2*c*d, d**2] ])
      para.append(twissTransform*paraInitial[i])

    return dct([('BETX',para[0].item(0)),
                ('ALFX',para[0].item(1)),
                ('GAMX',para[0].item(2)),
                ('BETY',para[1].item(0)),
                ('ALFY',para[1].item(1)),
                ('GAMY',para[1].item(2))])

  ## DISPERSION at location s in element nE
  # nE = number of element of interest, s = location of interest within element nE
  def getDisp(self, nE, s):
    # Get initial dispersion values DX, DPX, DY, DPY
    # If nE is first element, use start marker as initial; if nE any other element, use previous element
    if nE != 0:
      prevE = self.elems[nE-1]
    else:
      prevE = self.markers[0]

    if s == 0:
      return  dct([ ('DX',prevE.DX),
                    ('DPX',prevE.DPX),
                    ('DY',prevE.DY),
                    ('DPY',prevE.DPY) ])

    # Use same set-up as for new positions/angles to enable use of transport matrix as is
    # disp0 = [x=DX, px=DPX, y=DY, py=DPY, D=1, S=0]
    disp0 = mtrx( [ [prevE.DX],
                    [prevE.DPX],
                    [prevE.DY],
                    [prevE.DPY],
                    [1],
                    [0] ])

    # Copy element of interest and change length to location of interest, s. Get transport matrix assuming D=0
    # If nE is not DRIFT, QUADRUPOLE or DIPOLE, change element to DRIFT and recalculate transport matrix
    e = dct(self.elems[nE])
    e['L'] = s
    m = matrixForElement(e, 6)   # NOTE: Take order 6
    if m == None:
      e['KEYWORD']="DRIFT"
      m = matrixForElement(e, 6)
    m = m(d=0)

    # Multiply initial dispersions by transport matrix.
    # Equivalent to [n, n', 1]_new = [ [C S D], [C', S', D'], [0 0 1] ] * [n, n', 1]_old
    disp = m*disp0

    return dct([ ('DX',disp.item(0)),
                 ('DPX',disp.item(1)),
                 ('DY',disp.item(2)),
                 ('DPY',disp.item(3)) ])

  ## FIND ELEMENT: For a given s along the whole beamline, find in which element, i,  it is
  def findElem(self, s):
    for i in range(len(self.elems)):
      if s <= self.elems[i].S:
        return i

  ## NATURAL CHROMATICITY: Given a design Beta (or B*) and optionally
  ## the initial beta it calculates the natural chromaticity
  def getNatChrom(self, BetStar, BetY0 = None):
    newT = copy(self)
    newT.elems = []
    # Strip higher order elements from current twiss and store in a
    # new one
    for e in self.elems:
      if e.KEYWORD in ['DRIFT', 'QUADRUPOLE', 'SBEND'] and e.L != 0:
        newT.elems.append(e)
    # Calculate the map of the new twiss
    m = mapclass.Map2(newT)
    # For Xy, 001010 and Xy, 000110
    # Fr = gamma( (1+j+j') / 2 ) * gamma(...) * ... * gamma(3/2)
    Fr = 4.373354581906215 # gamma(1./2)**3 * gamma(3./2)**2
    # C = 2**((2+j+k+l+m+j'+k'+l'+m')/2) / pi**2.5
    C = 0.4573148512298903 # 8*pow(pi,-2.5)
    if BetY0 is None: BetY0 = self.markers[0].BETY
    return Fr*C*(m['y'][(0,0,1,0,1,0)]**2*BetY0/BetStar +
                 m['y'][(0,0,0,1,1,0)]**2/(BetY0*BetStar)).real

  ## OIDE EFFECT
  def oide(self,n=100):
    re = 2.817940325e-15
    lame = 3.861592678e-13
    emi = 0.2e-7
    betas = 17.92472388e-6
    gamma = 1500/0.000511
    coeff = 110*re*lame*(gamma**5)/(3*math.sqrt(6*math.pi))
    print coeff
    Kq = 0.1158931845985401
    Lq = 2.74
    ls = 3.5
    c = math.sqrt(Kq)*Lq
    def f(x):
      return math.sin(x) + math.sqrt(Kq)*ls*math.cos(x)
    def g(y):
      return (abs(math.sin(y) + math.sqrt(Kq)*ls*math.cos(y))**3)*(simpson(f,0,y,n))**2
    integral = simpson(g,0,c,n)
    ta = coeff*integral*(emi/(betas*gamma))**2.5
    tb = emi*betas/gamma
    ans = math.sqrt(ta+tb)
    return integral, ans


#########################
## Twiss functionality
#########################
def matrixForElement(e,order):
  try:
    r = None
    if e.KEYWORD == "DRIFT":
      r = DRIFT(order=order,**e)
    if e.KEYWORD == "QUADRUPOLE":
      if e.L != 0:
        if e.K1L > 0:
          r = QF(order=order,**e)
        else:
          r = QD(order=order,**e)
    if e.KEYWORD == "SBEND":
      r = DI(order=order,**e)
    return r
  except Exception as e:
    print "The Twiss object doesn't have the desired structure"
    print e
    exit()


def mapForElement(e,order):
  try:
    m = None
    if e.KEYWORD == "QUADRUPOLE":
      if e.L == 0:
        m = MUL(order=order,**e)
    if e.KEYWORD == "MULTIPOLE":
      if e.L == 0:
        m = MUL(order=order,**e)
    if e.KEYWORD in ["SEXTUPOLE", "OCTUPOLE", "DECAPOLE"]:
      if e.L == 0:
        m = MUL(order=order,**e)
      else:
        m = MULTHICK(order=order,**e)
    return m
  except Exception as e:
    print "The Twiss object doesn't have the desired structure"
    print e
    exit()
