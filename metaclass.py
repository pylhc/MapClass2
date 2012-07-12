from collections import namedtuple
from numpy import identity

from transport import *
from pytpsa import polmap

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

  #Calculate BETA, ALPHA & GAMMA at location s in element nE
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

    paraX0 = mtrx([ [betX0],
                    [alfX0],
                    [gamX0] ])

    paraY0 = mtrx([ [betY0],
                    [alfY0],
                    [gamY0] ])

    paraInitial = [paraX0, paraY0]

    # Copy element of interest and change its length to location of interest, s
    # Calculate transport matrix for this element assuming D=0
    e = self.elems[nE]                
    e['L'] = s                       
    eTransport = matrixForElement(e, 6)    # NOTE: Take order 6
    # If nE is not DRIFT, QUADRUPOLE or DIPOLE, change element to DRIFT and recalculate transport matrix
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
   
    return dict([('BETX',para[0].item(0)),
                 ('ALFX',para[0].item(1)),
                 ('GAMX',para[0].item(2)),
                 ('BETY',para[1].item(0)),
                 ('ALFY',para[1].item(1)),
                 ('GAMY',para[1].item(2))])

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
    if e.KEYWORD == "SEXTUPOLE" or \
       e.KEYWORD == "OCTUPOLE" or \
       e.KEYWORD == "DECAPOLE":
      if e.L == 0:
        m = MUL(order=order,**e)
      else:
        m = MULTHICK(order=order,**e)
    return m
  except Exception as e:
    print "The Twiss object doesn't have the desired structure"
    print e
    exit()
