from types import FunctionType
from numpy import matrix, nditer

from pytpsa import *

#########################
class mtrx(matrix):
#########################
  def __call__(self, *args, **kwargs):
    m = self.copy()
    for p in nditer(m, flags=['refs_ok'], op_flags=['readwrite']):
      if type(p.item()) is FunctionType:
        p[...] = p.item()(**kwargs)
    return m


## Polynomes declaration
X  = pol('x')
PX = pol('px')
Y  = pol('y')
PY = pol('py')
D  = pol('d')
S  = pol('s')

# Matrix with the incognites

U = matrix([ [X],
             [PX],
             [Y],
             [PY],
             [D],
             [S] ])

########################
# Transport matrices
########################

# DRIFT

def L(L,**args): return L/(1+D)

DRIFT = mtrx([ [1, L, 0, 0, 0, 0],
               [0, 1, 0, 0, 0, 0],
               [0, 0, 1, L, 0, 0],
               [0, 0, 0, 1, 0, 0],
               [0, 0, 0, 0, 1, 0],
               [0, 0, 0, 0, 0, 1] ])

# QUADRUPOLES

def Q11(L,K1L,**args): K = (K1L/L)/(1+D); return cos(L*sqrt(K))
def Q12(L,K1L,**args): K = (K1L/L)/(1+D); return (1/sqrt(K))*sin(L*sqrt(K))/(1+D)
def Q21(L,K1L,**args): K = (K1L/L)/(1+D); return -sqrt(K)*sin(L*sqrt(K))*(1+D)
def Q33(L,K1L,**args): K = (K1L/L)/(1+D); return cosh(L*sqrt(K))
def Q34(L,K1L,**args): K = (K1L/L)/(1+D); return (1/sqrt(K))*sinh(L*sqrt(K))/(1+D)
def Q43(L,K1L,**args): K = (K1L/L)/(1+D); return sqrt(K)*sinh(L*sqrt(K))*(1+D)

QF = mtrx([ [Q11, Q12, 0, 0, 0, 0],
            [Q21, Q11, 0, 0, 0, 0],
            [0, 0, Q33, Q34, 0, 0],
            [0, 0, Q43, Q33, 0, 0],
            [0, 0, 0, 0, 1, 0],
            [0, 0, 0, 0, 0, 1] ])

QD = mtrx([ [Q33, Q34, 0, 0, 0, 0],
            [Q43, Q33, 0, 0, 0, 0],
            [0, 0, Q11, Q12, 0, 0],
            [0, 0, Q21, Q11, 0, 0],
            [0, 0, 0, 0, 1, 0],
            [0, 0, 0, 0, 0, 1] ])
