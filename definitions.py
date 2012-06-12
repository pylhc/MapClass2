from pytpsa import polmap, mkpol
from numpy import matrix

XYZD=['x', 'px', 'y', 'py', 'd','s']

## Polynomes declaration
X,PX,Y,PY,D,S  = mkpol(str(XYZD).translate(None,"'[] "))

# Matrix with the incognites
def generateDefaultMatrix(order):
  X.order = PX.order = Y.order = PY.order = PY.order = D.order = S.order = order
  return matrix([ [X],
                  [PX],
                  [Y],
                  [PY],
                  [D],
                  [S] ])

def generateDefaultMap(order):
  return polmap(zip(XYZD,XYZD),order=order)
