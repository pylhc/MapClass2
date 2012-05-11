from pytpsa import pol
from numpy import matrix, nditer

#########################                                                                                
class mtrx(matrix):
#########################                                                                                
  def __call__(self, *args, **kwargs):
    for p in nditer(self, flags=['refs_ok'], op_flags=['readwrite']):
      if isinstance(p.item(), pol):
        p[...] = p.item()(**kwargs)



## Polynomes declaration
X  = pol('x')
PX = pol('px')
Y  = pol('y')
PY = pol('py')
D  = pol('d')
S  = pol('s')

###

L = pol('L')

# Matrix with the incognites

U = matrix([ [X],
             [PX],
             [Y],
             [PY],
             [D],
             [S] ])

# Transport matrices

DRIFT = mtrx([ [1, L, 0, 0, 0, 0],
               [0, 1, 0, 0, 0, 0],
               [0, 0, 1, L, 0, 0],
               [0, 0, 0, 1, 0, 0],
               [0, 0, 0, 0, 1, 0],
               [0, 0, 0, 0, 0, 1] ])
