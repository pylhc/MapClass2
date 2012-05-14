from pytpsa import *
from numpy import matrix, nditer

#########################                                                                                
class mtrx(matrix):
#########################                                                                                
  def __call__(self, *args, **kwargs):
    m = self.copy()
    for p in nditer(m, flags=['refs_ok'], op_flags=['readwrite']):
      if isinstance(p.item(), pol):
        p[...] = p.item()(**kwargs)
    return m


## Polynomes declaration
X  = pol('x')
PX = pol('px')
Y  = pol('y')
PY = pol('py')
D  = pol('d')
S  = pol('s')

###

L = pol('L')
K = pol('K1/(1+D)')

# Matrix with the incognites

U = matrix([ [X],
             [PX],
             [Y],
             [PY],
             [D],
             [S] ])

# Transport matrices

# PX = L/(1+D)
DRIFT = mtrx([ [1, L, 0, 0, 0, 0],
               [0, 1, 0, 0, 0, 0],
               [0, 0, 1, L, 0, 0],
               [0, 0, 0, 1, 0, 0],
               [0, 0, 0, 0, 1, 0],
               [0, 0, 0, 0, 0, 1] ])


QUADRUPOLE = mtrx([ [cos(L*sqrt(K)), (1/sqrt(K))*sin(L*sqrt(K)), 0, 0, 0, 0],
                    [-sqrt(K)*sin(L*sqrt(K)), cos(L*sqrt(K)), 0, 0, 0, 0],
                    [0, 0, cosh(L*sqrt(K)), (1/sqrt(K))*sinh(L*sqrt(K)), 0, 0],
                    [0, 0, sqrt(K)*sinh(L*sqrt(K)), cosh(L*sqrt(K)), 0, 0],
                    [0, 0, 0, 0, 1, 0],
                    [0, 0, 0, 0, 0, 1] ])
