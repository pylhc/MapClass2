from pytpsa import polmap, mkpol
from numpy import matrix

XYZD = ['x', 'px', 'y', 'py', 'd', 's']

## Polynomes declaration
X, PX, Y, PY, D, S = mkpol(str(XYZD).translate(None, "'[] "))


# Matrix with the incognites
def generateDefaultMatrix(order):
  X.order = PX.order = Y.order = PY.order = PY.order = D.order = S.order = order
  return matrix([[X],
                 [PX],
                 [Y],
                 [PY],
                 [D],
                 [S]])


def generateDefaultMap(order):
  return polmap(zip(XYZD, XYZD), order=order)

#########################
# Auxiliar elements
#########################

class dct(dict):
  '''
  Ladi da
  '''

  def __getattr__(self, attr):
    '''
    diii
    '''
    try:
      return self[attr]
    except:
      raise AttributeError("%r object has no attribute %r" % (type(self).__name__, attr))

  def __setattr__(self, attr, val):
    '''
    ddd
    '''
    self[attr] = val
