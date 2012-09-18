from numpy import identity
from copy import copy, deepcopy

from definitions import dct
from transport import *
import mapclass
from pytpsa import polmap
import math, re
from integration import simpson


#########################
class twiss2(dct):
#########################
  "Twiss parameters from madx output (with free choice of select items)"

  def __init__(self, filename='twiss'):
    self.elems = []
    # Markers (only takes 1st and last)
    self.markers = []

    f = open(filename, 'r')

    for line in f:
      # if ("@ " in line and "%le" in line) :    # FIX to take DPP %s
      if "@ " not in line and "@" in line:
        line = replace(line, "@", "@ ")

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

      if "$ " in line or "$\t" in line:
        types = splt[1:]

      if "@" not in line and "*" not in line and "%" not in line:
        vals = []
        for j in range(0, len(labels)):
          if "%hd" in types[j]:
            vals.append(int(splt[j]))
          if "%le" in types[j]:
            vals.append(float(splt[j]))
          if "s" in types[j]:
            vals.append(splt[j].strip('"'))
        e = dct(zip(labels, vals))
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

  def getBeta(self, nE, s):
    """
    Calculates beta, alpha, gamma at location s in element nE.

    :param int nE: number of element of interest - where the first element is 0
    :param float s: location of interest within element nE

    :return: dct of BETX, BETY, ALFX, ALFY, GAMX and GAMY
    """

    # Gets initial beta, alpha, gamma from previous element or start
    # marker (if nE=0)
    if nE != 0:
      prevE = self.elems[nE-1]
    else:
      prevE = self.markers[0]

    betX0 = prevE.BETX
    alfX0 = prevE.ALFX
    gamX0 = (1 + alfX0**2) / betX0

    betY0 = prevE.BETY
    alfY0 = prevE.ALFY
    gamY0 = (1 + alfY0**2) / betY0

    if s == 0:
      return  dct([('BETX', betX0),
                   ('ALFX', alfX0),
                   ('GAMX', gamX0),
                   ('BETY', betY0),
                   ('ALFY', alfY0),
                   ('GAMY', gamY0)])

    paraX0 = mtrx([[betX0],
                   [alfX0],
                   [gamX0]])

    paraY0 = mtrx([[betY0],
                   [alfY0],
                   [gamY0]])

    paraInitial = [paraX0, paraY0]

    # Copy element nE and change length to location of interest, s
    # Calculate transport matrix for this element assuming D=0
    # If nE is not DRIFT, QUADRUPOLE or DIPOLE, change element to
    # DRIFT and recalculate transport matrix
    e = dct(self.elems[nE])
    if e.L != 0:
      e.K1L = (e.K1L / e.L) * s
      e.ANGLE = (e.ANGLE / e.L) * s
    e.L = s
    eTransport = matrixForElement(e, 6)  # NOTE: Takes order 6
    if eTransport == None:
      e.KEYWORD ="DRIFT"
      eTransport = matrixForElement(e, 6)
    eTransport = eTransport(d=0)

    # Extract required elements from transport matrix to form transform matrix
    # i=0 for x-direction; i=1 for y-direction
    # a=C, b=S, c=C', d=S' according to standard notation
    para = []
    for i in [0, 1]:
      j = 2 * i
      a = eTransport.item((j, j))
      b = eTransport.item((j, j+1))
      c = eTransport.item((j+1, j))
      d = eTransport.item((j+1, j+1))
      twissTransform = mtrx([[a**2, -2*a*b, b**2],
                             [-a*c, b*c + a*d, -b*d],
                             [c**2, -2*c*d, d**2]])
      para.append(twissTransform*paraInitial[i])  # Calculate final values

    return dct([('BETX', para[0].item(0)),
                ('ALFX', para[0].item(1)),
                ('GAMX', para[0].item(2)),
                ('BETY', para[1].item(0)),
                ('ALFY', para[1].item(1)),
                ('GAMY', para[1].item(2))])

  def getDisp(self, nE, s):
    """
    Calculates dispersion at location s in element nE.

    :param int nE: number of element of interest
    :param float s: location of interest within element nE

    :return: dct of DX, DPX, DY and DPY
    """

    # Get initial dispersion values DX, DPX, DY, DPY from previous
    # element or start marker (if nE=0)
    if nE != 0:
      prevE = self.elems[nE-1]
    else:
      prevE = self.markers[0]

    if s == 0:
      return  dct([('DX', prevE.DX),
                   ('DPX', prevE.DPX),
                   ('DY', prevE.DY),
                   ('DPY', prevE.DPY)])

    # Set up initial "dispersion vector" such that multiplication by
    # transport matrix gives final dispersion function
    disp0 = mtrx([[prevE.DX],
                  [prevE.DPX],
                  [prevE.DY],
                  [prevE.DPY],
                  [1],
                  [0]])

    # Copy element nE and change length to location of interest, s
    # Calculate transport matrix for this element assuming D=0
    # If nE is not DRIFT, QUADRUPOLE or DIPOLE, change element to
    # DRIFT and recalculate transport matrix
    e = dct(self.elems[nE])
    if e.L != 0:
      e.K1L = (e.K1L / e.L) * s
      e.ANGLE = (e.ANGLE / e.L) * s
    e.L = s
    m = matrixForElement(e, 6)   # NOTE: Take order 6
    if m == None:
      e.KEYWORD = "DRIFT"
      m = matrixForElement(e, 6)
    m = m(d=0)

    # Calculate final values
    disp = m * disp0

    return dct([('DX', disp.item(0)),
                ('DPX', disp.item(1)),
                ('DY', disp.item(2)),
                ('DPY', disp.item(3))])

  def getPhase(self,nE,s):
    """
    Calculates phase at location s in element nE.

    :param int nE: number of element of interest - where the first element is 0
    :param float s: location of interest within element nE

    :return: dct of MUX and MUY
    """

    # Get initial phase values MUX and MUY from previous element
    # or start marker (if nE=0)
    if nE != 0:
      prevE = self.elems[nE-1]
    else:
      prevE = self.markers[0]

    if s == 0:
      return  dct([('MUX', prevE.MUX),
                   ('MUY', prevE.MUY)])

    para = self.getBeta(nE,s)

    # Copy element nE and change length to location of interest, s
    # Calculate transport matrix for this element assuming D=0
    # If nE is not DRIFT, QUADRUPOLE or DIPOLE, change element to
    # DRIFT and recalculate transport matrix
    e = dct(self.elems[nE])
    if e.L != 0:
      e.K1L = (e.K1L / e.L) * s
      e.ANGLE = (e.ANGLE / e.L) * s
    e.L = s
    m = matrixForElement(e, 6)  # NOTE: Takes order 6
    if m == None:
      e.KEYWORD = "DRIFT"
      m = matrixForElement(e, 6)
    m = m(d=0)

    # Calculate cos(delta Phi) and sin(delta Phi) in x and y planes
    xy = m.item((0,1)) / math.sqrt(prevE.BETX * para.BETX)
    xx = (math.sqrt(prevE.BETX) * m.item((0,0)) / math.sqrt(para.BETX)) - (prevE.ALFX * xy)

    yy = m.item((2,3)) / math.sqrt(prevE.BETY * para.BETY)
    yx = (math.sqrt(prevE.BETY) * m.item((2,2)) / math.sqrt(para.BETY)) - (prevE.ALFY * yy)

    thetaX = math.atan2(xy,xx)
    thetaY = math.atan2(yy, yx)
    if thetaX < 0:
      thetaX =+ 2 * math.pi
    if thetaY < 0:
      thetaY =+ 2 * math.pi

    return  dct([('MUX', thetaX / (2 * math.pi) + prevE.MUX),
                 ('MUY', thetaY / (2 * math.pi) + prevE.MUY)])

  def findElem(self, s):
    """
    Finds in which element a given location along the beamline is present

    :param float s: location along the whole beamline

    :return: integer number of element
    """

    for i in range(len(self.elems)):
      if s <= self.elems[i].S:
        return i

  def getNatChrom(self, BetStarX=None, BetX0=None, BetStarY=None, BetY0=None):
    """
    Returns the natural chromaticity of the beamline

    :param float BetStarX: design beta in x-direction
    :param float BetX0: initial beta in x-direction
    :param float BetStarY: design beta in y-direction
    :param float BetY0: initial beta in y-direction

    These parameters are optional and are otherwise read from the twiss markers
    """

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
    Fr = 4.373354581906215  # gamma(1./2)*3 * gamma(3./2)**2
    C = 0.4573148512298903  # 8*pow(pi,-2.5)

    if BetStarX is None: BetStarX = self.markers[1].BETX
    if BetX0 is None: BetX0 = self.markers[0].BETX
    if BetStarY is None: BetStarY = self.markers[1].BETY
    if BetY0 is None: BetY0 = self.markers[0].BETY
    #CHECK FOR X!!! JUST GUESSING
    return dct([('NChromX', Fr*C*(m['x'][(1,0,0,0,1,0)]**2 * BetX0 / BetStarX +
                                 m['x'][(0,1,0,0,1,0)]**2 / (BetX0 * BetStarX)).real),
                ('NChromY', Fr*C*(m['y'][(0,0,1,0,1,0)]**2 * BetY0 / BetStarY +
                                 m['y'][(0,0,0,1,1,0)]**2 / (BetY0 * BetStarY)).real)])

  def getChrom(self, s=None, s0=0, n=100):
    """
    Calculates chromaticity using -1/4pi * integral (beta*K) ds

    :param float s: end location along beamline
    :param float s0: start location along beamline (optional)
    :param int n: number of intervals for integration (optional)

    :returns: chromaticity between s0 and s
    """

    if s is None: s = self.markers[1].S

  ## CHECK: positive/negative signs on K for focus vs. defocus...
  ## Is this natural chromaticity also because it only considers quadrupoles?
  ## What about multipole quadrupoles?
    def fX(s):
      nE = self.findElem(s)
      e = self.elems[nE]
      ss = s - (e.S - e.L)
      bet = self.getBeta(nE, ss)
      if e.K1L != 0:
        return bet.BETX * -e.K1L / e.L  # Correct to make negative?
      else:
        return 0

    def fY(s):
      nE = self.findElem(s)
      e = self.elems[nE]
      ss = s - (e.S - e.L)
      bet = self.getBeta(nE, ss)
      if e.K1L != 0:
        return bet.BETY * e.K1L / e.L
      else:
        return 0

    return dct([('ChromX', -simpson(fX, s0, s, n) / (4 * math.pi)),
                ('ChromY', -simpson(fY, s0, s, n) / (4 * math.pi))])

  def oide(self, emi=2e-8, gamma=2.9354207436399e-6, n=100):
    """
    Returns delta(sigma^2) due to Oide Effect

    :param float emi: emittance
    :param float gamma: Lorentz factor = E/Eo = 1500/0.000511 for CLIC
    :param int n: number of intervals for integration (optional)
    """

    re = 2.817940325e-15
    lame = 3.861592678e-13
    betas = self.markers[1].BETY # Reads 6.77249e-5 from FFS
    # (betas = 17.92472388e-6 from mathematica nb)
    coeff = 110 * re * lame * gamma**5 / (3 * math.sqrt(6 * math.pi))

    # Read twiss object in reverse to find first DRIFT and QUADRUPOLE
    # to get ls, Lq and Kq
    for e in reversed(self.elems):
      if e.KEYWORD == 'DRIFT':
        ls = e.L
        break
    for e in reversed(self.elems):
      if e.KEYWORD == 'QUADRUPOLE':
        # Multiplied by 2 because final quadrupoles split in file
        Lq = 2 * e.L
        Kq = abs(e.K1L / e.L)
        break

    c = math.sqrt(Kq) * ls
    b = math.sqrt(Kq) * Lq

    # Define functions for integration
    def f(x):
      return math.sin(x) + c * math.cos(x)

    def g(y):
      return (abs(math.sin(y) + c * math.cos(y))**3) * simpson(f, 0, y, n)**2

    integral = simpson(g, 0, b, n)

    return coeff * integral * (emi / (betas * gamma))**2.5

  def getH(self, nE, s):
    """
    Returns H(s) function at location s in element nE

    :param int nE: number of element of interest
    :param float s: location of interest within element nE
    """

    para = self.getBeta(nE, s)
    disp = self.getDisp(nE, s)

    HX = (para.GAMX * disp.DX**2) + (2 * para.ALFX * disp.DX * disp.DPX) + (para.BETX * disp.DPX**2)
    HY = (para.GAMY * disp.DY**2) + (2 * para.ALFY * disp.DY * disp.DPY) + (para.BETY * disp.DPY**2)
    return dct([('HX', HX),
                ('HY', HY)])

  def sigmaBends(self, E, s=None, s0=0, n=100):
    """
    Returns delta(sigma^2) due to bends (dipoles)

    :param float E: energy
    :param float s: location of interest along beamline (optional)
    :param float s0: start location along beamline (optional)
    :param int n: number of intervals for integrations (optional)
    """

    if s is None:
      s = self.markers[1].S
      endPhase = self.markers[1].MUX
    else:
      nE = self.findElem(s)
      e = self.elems[nE]
      ss = s - (e.S - e.L)
      endPhase = self.getPhase(nE, ss)

    # Calculates H*G^3 cos(phi)^2 at location s along the beamline
    def f(s):
      nE = self.findElem(s)
      e = self.elems[nE]

      # Calculate function only if element is dipole (i.e. ANGLE not 0)
      if e.ANGLE != 0:
        # ss from beginning of element nE == s from beginning of beamline
        ss = s - (e.S - e.L)
        para = self.getBeta(nE,ss)
        disp = self.getDisp(nE,ss)
        alpha = math.atan(-para.ALFX - para.BETX * disp.DPX / disp.DX)
        Phi = (endPhase - self.getPhase(nE,ss).MUX) * 2 * math.pi + alpha
        cosPhi = (math.cos(Phi) ** 2).real
        H = self.getH(nE, ss)
        P = abs(e.L / e.ANGLE)
        return H.HX * cosPhi / P**3
      else:
        return 0

    c2 = 4.13e-11  # m^2(GeV)^-5
    coeff = c2 * E**5 * self.markers[1].BETX
    return coeff * simpson(f, s0, s, n)

  def sigmaBends2(self, E, s=None, s0=0, n=20):
    """
    Returns delta(sigma^2) due to bends (dipoles)

    :param float E: energy
    :param float s: location of interest along beamline (optional)
    :param float s0: start location along beamline (optional)
    :param int n: number of intervals for integrations (optional)
    """

    if s is None:
      s = self.markers[1].S
      nELast = len(self.elems)-1
      endPhase = self.markers[1].MUX
      ss = self.elems[-1].L
    else:
      nELast = self.findElem(s)
      last = self.elems[nELast]
      ss = s - (last.S - last.L)
      endPhase = self.getPhase(nELast, ss).MUX

    if s0 != 0:
      nEFirst = self.findElem(s0)
      first = self.elems[nEFirst]
      ss0 = s0 - (first.S - first.L)
      rangeElems = xrange(nEFirst, len(self.elems))
    else:
      nEFirst = 0
      rangeElems = xrange(len(self.elems))

    # Calculates H*G^3 cos(phi)^2 at location s along the beamline
    def wrap(nE):
      def f(ss):
        para = self.getBeta(nE,ss)
        disp = self.getDisp(nE,ss)
        if disp.DX == 0: return 0
        alpha = math.atan(-para.ALFX - para.BETX * disp.DPX / disp.DX)
        Phi = (endPhase - self.getPhase(nE,ss).MUX) * 2 * math.pi + alpha
        cosPhi = (math.cos(Phi) ** 2).real
        H = self.getH(nE, ss)
        P = abs(e.L / e.ANGLE)
        return H.HX * cosPhi / P**3
      return f

    c2 = 4.13e-11  # m^2(GeV)^-5
    coeff = c2 * E**5 * self.markers[1].BETX

    total = 0
    for i in rangeElems:
      e = self.elems[i]
      if e.KEYWORD == 'SBEND':
        if i == nEFirst:
          total +=  coeff * simpson(wrap(i), ss0, e.L, n)
        elif i == nELast:
          total += coeff * simpson(wrap(i), 0, ss, n)
        else:
          total += coeff * simpson(wrap(i), 0, e.L, n)
      if i == nELast: return total

  def stripLine(self):
    """
    Returns a new twiss object with the monitors and markers and matrices removed
    """

    t = deepcopy(self)
    t.elems = [e for e in t.elems if e.KEYWORD not in ["MARKER", "MONITOR", "MATRIX"]]
    return t

  def mergeElems(self):
    """
    Returns a new twiss object with adjacent elements combined if they have the
    same KEYWORD, bending radius and KnL.
    """

    t = deepcopy(self)
    i = 0
    while i < len(t.elems) - 1:
      curr = t.elems[i]
      nxt = t.elems[i+1]
      # Make subdictionaries of KEYWORD, bending radius and all of the strength
      # parameters KnL for quick comparison
      currSub = dict((k, v) for k, v in curr.iteritems() if re.match("K\d+L", k) or k in ["KEYWORD"])
      nxtSub =  dict((k, v) for k, v in nxt.iteritems() if re.match("K\d+L", k) or k in ["KEYWORD"])
      if currSub["KEYWORD"] == "SBEND" and nxtSub["KEYWORD"] == "SBEND":
        currSub["RHO"] = curr.ANGLE / curr.L
        nxtSub["RHO"] = nxt.ANGLE / nxt.L
      # If subdictionaries are equal change KnL, ANGLE and L of the
      # second element and delete the first one.
      if currSub == nxtSub:
        if nxt.L != 0:
          for k,knl in nxt.iteritems():
            if re.match("K\d+L", k) and knl != 0:
              nxt[k] = (knl / nxt.L) * (curr.L + nxt.L)
          nxt.ANGLE = (nxt.ANGLE / nxt.L) * (curr.L + nxt.L)
          nxt.L = nxt.L + curr.L
        del t.elems[i]
      # i only increments when subdictionaries not equal because upon
      # deletion of an element len(t.elems) decreases by 1 for the loop
      else: i = i + 1
    return t

  def alterElem(self, nE, dL=0, dPos=0):
    """
    Returns a new twiss object having altered the length or position (or both)
    of an element.

    Changes the twiss parameters (BETX, BETY, ALFX, ALFY), the dispersion
    (DX, DPX, DY, DPY) the phase (MUX, MUY) and the strength parameters KnL
    accordingly. Other parameters will maintain the values from the original
    twiss object and, therefore, may be incorrect.

    :param int nE: number of element whose properties are to be altered
    :param float dL: change in length (e.g. dL = 2 adds 1 unit to each end of the element)
    :param float dPos: change in position (e.g. dPos = 2 moves it 2 units forward along the line)
    """

    # Tests that element is not first or last in the line
    if nE < 1 or nE > len(self.elems) - 2:
      print "Element out of bounds"
      return

    # Tests that element is surrounded by drifts
    prev = self.elems[nE-1]
    nxt = self.elems[nE+1]
    for k in [prev, nxt]:
      if k.KEYWORD != "DRIFT":
        print "Element not surrounded by drifts"
        return

    # Tests that dL not longer than surrounding drift space
    if dL > (prev.L + nxt.L):
      print "dL too long"
      return

    # Tests that dPos does not exceed available drift space
    if dPos < 0:
      if abs(dPos) > prev.L:
        print "dPos out of range"
        return
    if dPos > 0:
      if dPos > nxt.L:
        print "dPos out of range"
        return

    # Makes a copy of the twiss object
    t = deepcopy(self)
    prev = t.elems[nE-1]
    curr = t.elems[nE]
    nxt = t.elems[nE+1]

    # Modifies L and S of twiss copy according to length change, dL
    # Length is added/subtracted symmetrically from each side of the element
    prev.L = prev.L - dL / 2.0
    if curr.L != 0:
      for k,knl in curr.iteritems():
        if re.match("K\d+L", k) and knl != 0:
          curr[k] = (knl / curr.L) * (dL + curr.L)
      curr.ANGLE = (curr.ANGLE / curr.L) * (dL + curr.L)
    curr.L = curr.L + dL
    nxt.L = nxt.L - dL / 2.0

    prev.S = prev.S - dL / 2.0
    curr.S = curr.S + dL / 2.0

    # Modifies L and S of twiss copy according to position change, dPos
    prev.L = prev.L + dPos
    nxt.L = nxt.L - dPos

    prev.S = prev.S + dPos
    curr.S = curr.S + dPos

    # If Quadrupole or Dipole, change twiss parameters, dispersion and
    # phase till the end of the line
    # If any other element, change till the end of the second drift
    # (as they are modelled as drifts anyway)
    if curr.KEYWORD in ["QUADRUPOLE", "SBEND"]:
      end = len(t.elems)
    else:
      end = nE + 1

    for i in range(nE-1, end):
      e = t.elems[i]
      para = t.getBeta(i,e.L)
      disp = t.getDisp(i,e.L)
      mu = t.getPhase(i,e.L)
      for k in para.keys():
        e[k] = para[k]
      for k in disp.keys():
        e[k] = disp[k]
      for k in mu.keys():
        e[k] = mu[k]

    return t


#########################
## Twiss functionality
#########################
def matrixForElement(e, order):
  try:
    r = None
    if e.KEYWORD == "DRIFT":
      r = DRIFT(order=order, **e)
    if e.KEYWORD == "QUADRUPOLE":
      if e.L != 0:
        if e.K1L > 0:
          r = QF(order=order, **e)
        else:
          r = QD(order=order, **e)
    if e.KEYWORD == "SBEND":
      r = DI(order=order, **e)
    return r
  except Exception as e:
    print "The Twiss object doesn't have the desired structure"
    print e
    exit()


def mapForElement(e, order):
  try:
    m = None
    if e.KEYWORD == "QUADRUPOLE":
      if e.L == 0:
        m = MUL(order=order, **e)
    if e.KEYWORD == "MULTIPOLE":
      if e.L == 0:
        m = MUL(order=order, **e)
    if e.KEYWORD in ["SEXTUPOLE", "OCTUPOLE", "DECAPOLE"]:
      if e.L == 0:
        m = MUL(order=order, **e)
      else:
        m = MULTHICK(order=order, **e)
    return m
  except Exception as e:
    print "The Twiss object doesn't have the desired structure"
    print e
    exit()
