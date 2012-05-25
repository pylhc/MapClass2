from pytpsa import polmap, mkpol
from numpy import matrix

XYZD=['x', 'px', 'y', 'py', 'd','s']

## Polynomes declaration
X,PX,Y,PY,D,S  = mkpol(str(XYZD).translate(None,"'[] "))

# Matrix with the incognites
U = matrix([ [X],
             [PX],
             [Y],
             [PY],
             [D],
             [S] ])

def generateDefaultMap():
  return polmap(**dict(zip(XYZD,XYZD)))
