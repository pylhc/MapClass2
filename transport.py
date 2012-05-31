from types import FunctionType
from numpy import matrix, nditer
from math import factorial

from definitions import *
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


#### Helper functions ####
def matrixToMap(m, vrs):
  mp = polmap()
  for i in range(0,len(vrs)):
    mp[vrs[i]] = m.item(i)
  return mp



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

def Q11(L,K1L,**args): K = abs(K1L/L)/(1+D); return cos(L*sqrt(K))
def Q12(L,K1L,**args): K = abs(K1L/L)/(1+D); return (1/sqrt(K))*sin(L*sqrt(K))/(1+D)
def Q21(L,K1L,**args): K = abs(K1L/L)/(1+D); return -sqrt(K)*sin(L*sqrt(K))*(1+D)
def Q33(L,K1L,**args): K = abs(K1L/L)/(1+D); return cosh(L*sqrt(K))
def Q34(L,K1L,**args): K = abs(K1L/L)/(1+D); return (1/sqrt(K))*sinh(L*sqrt(K))/(1+D)
def Q43(L,K1L,**args): K = abs(K1L/L)/(1+D); return sqrt(K)*sinh(L*sqrt(K))*(1+D)

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

# Different implementation for thin quadrupole
# Now the implementation for thin multipole is used
# def QT21(K1L,**args): return -K1L
# def QT43(K1L,**args): return K1L

# QTHIN = mtrx([ [1, 0, 0, 0, 0, 0],
#                [QT21, 1, 0, 0, 0, 0],
#                [0, 0, 1, 0, 0, 0],
#                [0, 0, QT43, 1, 0, 0],
#                [0, 0, 0, 0, 1, 0],
#                [0, 0, 0, 0, 0, 1] ])

# DIPOLES

def D11(L,ANGLE,**args): THETA = ANGLE/sqrt(1+D); return cos(THETA)
def D12(L,ANGLE,**args): THETA = ANGLE/sqrt(1+D); P = (L/ANGLE)/sqrt(1+D); return P*sin(THETA)
def D15(L,ANGLE,**args): THETA = ANGLE/sqrt(1+D); P = L/ANGLE; return P*(1-cos(THETA))
def D21(L,ANGLE,**args): THETA = ANGLE/sqrt(1+D); P = (L/ANGLE)/sqrt(1+D); return -(1/P)*sin(THETA)
def D25(L,ANGLE,**args): THETA = ANGLE/sqrt(1+D); return sin(THETA)*sqrt(1+D)
def D34(L,**args): return L/(1+D)

DI = mtrx([ [D11, D12, 0, 0, D15, 0],
            [D21, D11, 0, 0, D25, 0],
            [0, 0, 1, D34, 0, 0],
            [0, 0, 0, 1, 0, 0],
            [0, 0, 0, 0, 1, 0],
            [0, 0, 0, 0, 0, 1] ])

# MULTIPOLES

J = complex(0,1)

def EQ(n): return (1./factorial(n))*(X+J*Y)**n

def Bx(N,KnL): r,i = EQ(N).separateComplex(); return KnL*-r
def By(N,KnL): r,i = EQ(N).separateComplex(); return KnL*i

# It doesn't go any higher than 'decapole' but you can simply add it if you are
# going to use it
def MUL(K1L,K2L,K3L,K4L,**args):
  m = generateDefaultMap()
  m[XYZD[1]] = PX + Bx(1,K1L) + Bx(2,K2L) + Bx(3,K3L) + Bx(4,K4L)
  m[XYZD[3]] = PY + By(1,K1L) + By(2,K2L) + By(3,K3L) + By(4,K4L)
  return m


# THICK MULTIPOLES - Simpsons approximation
#
# Ci/k = (41/840, 9/35, 9/280, 34/105, 9/280, 9/35, 41/840)
#
# Note: The list is symetric

CiK = [0.04880952380952381, 0.2571428571428571, 0.03214285714285714, 0.3238095238095238,
       0.03214285714285714, 0.2571428571428571, 0.04880952380952381]

def MULTHICK(K1L,K2L,K3L,K4L,L,**args):
  Li = L/6 # Divide the original length in 6 splits
  m = generateDefaultMap()

  for ck in CiK[:-1]:
    m = MUL(K1L*ck,K2L*ck,K3L*ck,K4L*ck,**args) * m
    m = matrixToMap(DRIFT(L=Li,**args) * U, XYZD) * m

  m = MUL(K1L*CiK[-1],K2L*CiK[-1],K3L*CiK[-1],K4L*CiK[-1],**args) * m
  return m
