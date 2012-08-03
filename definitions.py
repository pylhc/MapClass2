from pytpsa import polmap, mkpol
from numpy import matrix


XYZD = ['x', 'px', 'y', 'py', 'd', 's']

## Polynomes declaration
X, PX, Y, PY, D, S = mkpol(str(XYZD).translate(None, "'[] "))


def generateDefaultMatrix(order):
  """
  A default matrix is a single column matrix with a
  simple polynomial per dimension.
  """
  X.order = PX.order = Y.order = PY.order = PY.order = D.order = S.order = order
  return matrix([[X],
                 [PX],
                 [Y],
                 [PY],
                 [D],
                 [S]])


def generateDefaultMap(order):
  """
  A default map is a map where all the polynomials are simple;
  this could be seen as an "identity" map, where if you combine (*) it with another
  map you obtain the previous one (i.e. defaultMap * mapA = mapA ).

  """
  return polmap(zip(XYZD, XYZD), order=order)

#########################
# Auxiliar elements
#########################

class dct(dict):
  '''
  Extends the dictionary class to allow the access of elements of
  the dictionary as attributes (i.e. d['key1'] vs. d.key1).
  '''

  def __getattr__(self, attr):
    '''
    Query the dictionary when an attribute is requested
    '''
    try:
      return self[attr]
    except:
      raise AttributeError("%r object has no attribute %r" % (type(self).__name__, attr))

  def __setattr__(self, attr, val):
    '''
    Set an element of the dictionary when an attribute is set
    '''
    self[attr] = val
