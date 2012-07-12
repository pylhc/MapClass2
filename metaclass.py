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

      if "@" not in line and "*" not in line and "$" not in line:
        vals = []
        for j in range(0,len(labels)):
          if "%hd" in types[j]:
            vals.append(int(splt[j]))
          if "%le" in types[j]:
            vals.append(float(splt[j]))
          if "s" in types[j]:
            vals.append(splt[j].strip('"'))
        self.elems.append(dct(zip(labels,vals)))

    f.close()
    try:
      labels
      types
    except:
      print "From Metaclass: Bad format or empty file ", filename
      print "Leaving Metaclass"
      exit()

  # Add a BEAMBEAM element - by default it adds at the end
  def addBeamBeam(self, name, K, sigmax, sigmay, xm, ym, position = None, t = "ROUNDBB"):
    e = dct({'KEYWORD': t, 'NAME':name, 'K':K, 'SIGMAX':sigmax, 'SIGMAY':sigmay, 'XM':xm, 'YM':ym})
    if position == None: position = len(self.elems)
    self.elems.insert(position, e)

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
    if e.KEYWORD == "ROUNDBB":
      m = ROUNDBB(order=order, **e)
    return m
  except Exception as e:
    print "The Twiss object doesn't have the desired structure"
    print e
    exit()
