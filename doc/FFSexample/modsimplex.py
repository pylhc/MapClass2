#!/usr/bin/env python
# 
# Copyright (c) 2001 Vivake Gupta (v@omniscia.org).  All rights reserved.
# 
# This program is free software; you can redistribute it and/or
# modify it under the terms of the GNU General Public License as
# published by the Free Software Foundation; either version 2 of the
# License, or (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program; if not, write to the Free Software
# Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307
# USA
#
# This software is maintained by Vivake (v@omniscia.org) and is available at:
#     http://www.omniscia.org/~vivake/python/Simplex.py

# Modified (debugged?) 7/16/2004 Michele Vallisneri (vallis@vallis.org)

""" Simplex - a regression method for arbitrary nonlinear function minimization

Simplex minimizes an arbitrary nonlinear function of N variables by the
Nelder-Mead Simplex method as described in:

Nelder, J.A. and Mead, R., "A Simplex Method for Function Minimization",
   Computer Journal, Vol. 7, 1965, pp. 308-313.

It makes no assumptions about the smoothness of the function being minimized.
It converges to a local minimum which may or may not be the global minimum
depending on the initial guess used as a starting point.
"""

#from Numeric import *
from string import split
#from MLab import mean
#from LinearAlgebra import singular_value_decomposition
#from LinearAlgebra import *
#import sys
#from random import gauss, seed
#from FFT import real_fft

import math
import copy
import sys
from os import system

from ctypes import *

cdll.LoadLibrary("../../libs/boost_1_53_0/libboost_python.so.1.53.0")
sys.path.append('../../')
sys.path.append('../../libs')
sys.path.append('../../doc/FFSexample')

class Simplex:
    def __init__(self, testfunc, guess, increments, kR = -1, kE = 2, kC = 0.5):
        """
        Initializes the simplex.
        INPUTS
        ------
        testfunc      the function to minimize
        guess[]       a list containing initial guesses
        increments[]  a list containing increments, perturbation size
        kR            step size/learning constant (sigma = ??)
        kE            expansion constant   (gamma = 2.0)
        kC            contraction constant (beta  = 0.5)
        """
        self.testfunc = testfunc
        self.guess = guess
        self.increments = increments
        self.kR = kR
        self.kE = kE
        self.kC = kC
        self.numvars = len(self.guess)
        self.simplex = []

        self.lowest = -1
        self.highest = -1
        self.secondhighest = -1

        self.errors = []
        self.currenterror = 0
        
        self.g = [0]*self.numvars
        self.f = open('simplexprog.dat', 'w')
        # Initialize vertices
        # MV: the first vertex is just the initial guess
        #     the other N vertices are the initial guess plus the individual increments
        #     the last two vertices will store the centroid and the reflected point
        #     the compute errors at the ... vertices
        
        print "Calculating simplex"
        for vertex in range(0, self.numvars + 3):
            self.simplex.append(copy.copy(self.guess))

        for vertex in range(0, self.numvars + 1):
            for x in range(0, self.numvars):
                if x == (vertex - 1):
                    self.simplex[vertex][x] = self.guess[x] + self.increments[x]
            self.errors.append(0)    
        print "Calculating function values at vertices"
        self.calculate_errors_at_vertices()
        
    def minimize(self, epsilon = 1e-20, maxiters = 1000, monitor = 1):
    
        from math import sqrt
        
        """Walks to the simplex down to a local minima.
        INPUTS
        ------
        epsilon       convergence requirement
        maxiters      maximum number of iterations
        monitor       if non-zero, progress info is output to stdout  

        OUTPUTS
        -------
        an array containing the final values
        lowest value of the error function
        number of iterations taken to get here
        """
        
        iter = 0

        for iter in range(0, maxiters):
        
            # Identify highest and lowest vertices
            #print "Identifying highest/lowest vertices"
            self.highest = 0
            self.lowest = 0
            for vertex in range(1, self.numvars + 1):
                if self.errors[vertex] > self.errors[self.highest]:
                    self.highest = vertex
                if self.errors[vertex] < self.errors[self.lowest]:
                    self.lowest = vertex
            print "highest: vertex ", self.highest, '\t', self.simplex[self.highest]
            print "lowest: vertex ", self.lowest, '\t', self.simplex[self.lowest]
                         
            # Move highest vertex to end
            print "Moving highest vertex to end"
            temp = self.guess
            if self.highest != self.numvars:
                temp = self.simplex[self.highest]
                if self.lowest == self.numvars:
                    self.lowest = self.highest
                    self.errors[self.lowest] = self.errors[self.numvars]
                self.simplex[self.highest] = self.simplex[self.numvars]
                self.simplex[self.numvars] = temp
                self.highest = self.numvars
            self.calculate_errors_at_vertices()
            print "points: ", self.simplex
            print "function values: ", self.errors
            
            # Identify second-highest vertex
            #print "Identifying second-highest vertex"
            self.secondhighest = self.lowest
            for vertex in range(0, self.numvars + 1):
                if vertex == self.highest:
                    continue
                elif vertex == self.secondhighest:
                    continue
                elif self.errors[vertex] > self.errors[self.secondhighest]:
                    self.secondhighest = vertex
            print "2nd highest: vertex ", self.secondhighest, self.simplex[self.secondhighest]

            # Test for convergence:
            #   compute the average merit figure (ugly)
            print "Testing for convergence"
            S = 0.0
            for vertex in range(0, self.numvars + 1):
                S = S + self.errors[vertex]
            F2 = S / (self.numvars + 1)

            #   compute the std deviation of the merit figures (ugly)

            S1 = 0.0
            for vertex in range(0, self.numvars + 1):
                S1 = S1 + (self.errors[vertex] - F2)**2
            T = math.sqrt(S1 / self.numvars)
            
            # Optionally, print progress information

            if monitor:
                print '\r' + 72 * ' ',
                print '\rIteration = %d   Best = %f   Worst = %f' % \
                      (iter,self.errors[self.lowest],self.errors[self.highest]),
                sys.stdout.flush()
                
            if T <= epsilon: #and self.testfunc(self.simplex[self.highest]) <= epsilon:
                # We converged!  Break out of loop!
                print "We converged! break"
                break;
            else:
                # Didn't converge.  Keep crunching.
                
                # Calculate xs
                # store xs in element N + 1
                print "\nCalculating xs"
                for x in range(0, self.numvars):
                    self.simplex[self.numvars + 1][x] = self.simplex[x][x]
                    if x % 2 == 0 and self.simplex[x][x] == self.simplex[x + 1][x]:
                        self.simplex[self.numvars + 1][x] = self.simplex[x + 1][x] + 1e-5
                    elif x % 2 == 1 and self.simplex[x][x] == self.simplex[x - 1][x]:
                        self.simplex[self.numvars + 1][x] = self.simplex[x - 1][x] + 1e-5
                print "xs = ", self.simplex[self.numvars + 1]

                self.gradient_simplex()

                # reflect the simplex across with pseudo-gradient
                # store reflected point in elem. N + 1 (and self.guess)
                
                self.reflect_simplex()
                self.currenterror = self.testfunc(self.guess)

                if self.currenterror < self.errors[self.highest]:
                    self.accept_reflected_point()

                if self.currenterror <= self.errors[self.lowest]:
                    self.expand_simplex()
                    self.currenterror = self.testfunc(self.guess)

                    # at this point we can assume that the highest
                    # value has already been replaced once
                    if self.currenterror < self.errors[self.highest]:
                        self.accept_expanded_point()
                elif self.currenterror >= self.errors[self.secondhighest]:
                    # worse than the second-highest, so look for
                    # intermediate lower point

                    self.contract_simplex()
                    self.currenterror = self.testfunc(self.guess)

                    if self.currenterror < self.errors[self.highest]:
                        self.accept_contracted_point()
                    else:
                        self.multiple_contract_simplex()
 
        # Either converged or reached the maximum number of iterations.
        # Return the lowest vertex and the currenterror.
        for vertex in range(1, self.numvars + 1):
            if self.errors[vertex] > self.errors[self.highest]:
                self.highest = vertex
            if self.errors[vertex] < self.errors[self.lowest]:
                self.lowest = vertex

        for x in range(0, self.numvars):
            self.guess[x] = self.simplex[self.lowest][x]
        self.currenterror = self.errors[self.lowest]
        return self.guess, self.currenterror, iter

    def gradient_simplex(self):
        # Calculate pseudo-gradients
        print "Calculating quasi-gradients"
        for x in range(0, self.numvars):
            if x % 2 == 1:
                self.g[x] = (self.testfunc(self.simplex[x-1]) - self.testfunc(self.simplex[self.numvars+1]))/(self.simplex[x-1][x] - self.simplex[self.numvars+1][x])
            else:
                 self.g[x] = (self.testfunc(self.simplex[x+1]) - self.testfunc(self.simplex[self.numvars+1]))/(self.simplex[x+1][x] - self.simplex[self.numvars+1][x])
        print "G = ", self.g
        return   

    # same as expand, but with alpha < 1; kC = 0.5 fine with NR
    def contract_simplex(self):
        print "contracting"
        for x in range(0, self.numvars):
            #self.guess[x] = self.kC * self.simplex[self.highest][x] + (1 - self.kC) * self.simplex[self.numvars+1][x]
            self.guess[x] = self.kC * self.simplex[self.numvars + 2][x] + (1 - self.kC) * self.simplex[self.lowest][x]
        print "Trying ", self.guess
        return

    # expand: if R is refelcted point and B is best point, alpha-expansion is
    #         (1 - alpha)*B + alpha*R; default alpha is 2.0; agrees with NR
    def expand_simplex(self):
        print "expanding"
        for x in range(0, self.numvars):
            #self.guess[x] = self.kE * self.guess[x] + (1 - self.kE) * self.simplex[self.numvars + 1][x]
            self.guess[x] = self.kE * self.simplex[self.numvars + 2][x] + (1 - self.kE) * self.simplex[self.lowest][x]
        print "Trying ", self.guess
        return

    # reflect: if G is pseudo-gradients and Q is best point, reflection is B - sigma * G,
    #          which is achieved for kR = -1 (default value); agrees with NR
    def reflect_simplex(self):
        # loop over variables
        print "reflecting"
        for x in range(0, self.numvars):
            #self.guess[x] = self.kR * self.simplex[self.highest][x] + (1 - self.kR) * self.simplex[self.numvars+1][x]
            self.guess[x] = -1*self.kR * self.simplex[self.lowest][x] -0.5 * (1 - self.kR) * self.g[x]
            self.simplex[self.numvars + 2][x] = self.guess[x]
        print "Trying ", self.guess            
        # store reflected point in elem. N + 2

             
        return

    # multiple contraction: around the lowest point; agrees with NR

    def multiple_contract_simplex(self):
        print "contracting simplex"
        for vertex in range(0, self.numvars + 1):
            if vertex == self.lowest:
                continue
            for x in range(0, self.numvars):
                if self.simplex[vertex][x] == 0.5 * (self.simplex[vertex][x] + self.simplex[self.lowest][x]):
                    self.simplex[vertex][x] = (self.simplex[vertex][x] + self.simplex[self.lowest][x])*(vertex + 1)
                else:
                    self.simplex[vertex][x] = 0.5 * (self.simplex[vertex][x] + self.simplex[self.lowest][x])
        self.calculate_errors_at_vertices()
        return

    def accept_contracted_point(self):
        print "accepting contracted point: ", self.guess
        for x in range(0, self.numvars):
            self.simplex[self.highest][x] = self.guess[x]
        self.errors[self.highest] = self.currenterror
        # print vertices, xs, and gradients            
        for vertex in range(0, self.numvars + 1):
            print "vertex ", vertex, " = ", self.simplex[vertex], "function = ", self.errors[vertex]
        print "xs = ",self.simplex[self.numvars+1], "G = ", self.g
        # write simplex to file 
            
        for vertex in range(0, self.numvars + 3):
            x = str(self.simplex[vertex][0])
            y = str(self.simplex[vertex][1])
            self.f.write(x)
            self.f.write('\t')
            self.f.write(y)
            self.f.write('\t')
        self.f.write(str(self.g[0]))
        self.f.write('\t')
        self.f.write(str(self.g[1]))
        self.f.write('\t')
        self.f.write(str(self.lowest))
        self.f.write('\n')
        return

    def accept_expanded_point(self):
        print "accepting expanded point: ", self.guess
        for x in range(0, self.numvars):
            self.simplex[self.highest][x] = self.guess[x]
        self.errors[self.highest] = self.currenterror
        # print vertices, xs, and gradients            
        for vertex in range(0, self.numvars + 1):
            print "vertex ", vertex, " = ", self.simplex[vertex], "function = ", self.errors[vertex]
        print "xs = ",self.simplex[self.numvars+1], "G = ", self.g
        # write simplex to file 
            
        for vertex in range(0, self.numvars + 3):
            x = str(self.simplex[vertex][0])
            y = str(self.simplex[vertex][1])
            self.f.write(x)
            self.f.write('\t')
            self.f.write(y)
            self.f.write('\t')
        self.f.write(str(self.g[0]))
        self.f.write('\t')
        self.f.write(str(self.g[1]))
        self.f.write('\t')
        self.f.write(str(self.lowest))
        self.f.write('\n')
        return

    def accept_reflected_point(self):
        print "accepting reflected point: ", self.guess
        for x in range(0, self.numvars):
            self.simplex[self.highest][x] = self.simplex[self.numvars + 2][x]
        self.errors[self.highest] = self.currenterror
        # print vertices, xs, and gradients            
        for vertex in range(0, self.numvars + 1):
            print "vertex ", vertex, " = ", self.simplex[vertex], "function = ", self.errors[vertex]
        print "xs = ",self.simplex[self.numvars+1], "G = ", self.g
        # write simplex to file 
            
        for vertex in range(0, self.numvars + 3):
            x = str(self.simplex[vertex][0])
            y = str(self.simplex[vertex][1])
            self.f.write(x)
            self.f.write('\t')
            self.f.write(y)
            self.f.write('\t')
        self.f.write(str(self.g[0]))
        self.f.write('\t')
        self.f.write(str(self.g[1]))
        self.f.write('\t')
        self.f.write(str(self.lowest))
        self.f.write('\n')
        return

    def calculate_errors_at_vertices(self):
        for vertex in range(0, self.numvars + 1):
            if vertex == self.lowest:
                # compute the error unless we're the lowest vertex
                continue
            for x in range(0, self.numvars):
                self.guess[x] = self.simplex[vertex][x]
            self.currenterror = self.testfunc(self.guess)
            self.errors[vertex] = self.currenterror
        # print vertices, xs, and gradients            
        for vertex in range(0, self.numvars + 1):
            print "vertex ", vertex, " = ", self.simplex[vertex], "function = ", self.errors[vertex]
        print "xs = ",self.simplex[self.numvars+1], "G = ", self.g
        # write simplex to file 
            
        for vertex in range(0, self.numvars + 3):
            x = str(self.simplex[vertex][0])
            y = str(self.simplex[vertex][1])
            self.f.write(x)
            self.f.write('\t')
            self.f.write(y)
            self.f.write('\t')
        self.f.write(str(self.g[0]))
        self.f.write('\t')
        self.f.write(str(self.g[1]))
        self.f.write('\t')
        self.f.write(str(self.lowest))
        self.f.write('\n')
        
        return

#########################
def sigma(deltafamilies):
#########################
    global betx, bety, gamma, ex, ey, SF1, SD0
    #not necessary since no longer calling function
    
    from metaclass2 import twiss2
    from mapclass import Map2
    from math import sqrt
    
    order = 3
    
    gaussian = False    
    betx=64.99988501
    bety=17.99971417
    gamma=3e6
    ex=68e-8
    ey=2e-8
    
    t0 = twiss2('assets/ffs.twiss')
    t0 = t0.stripLine()
    t0 = t0.mergeElems()
 
    #print 'Deltafamilies = ', deltafamilies
    t0.elems[24].K2L = deltafamilies[0]
    t0.elems[30].K2L = deltafamilies[1]
    t0.elems[31].K2L = deltafamilies[2]
    t0.elems[32].K2L = deltafamilies[3]
    t0.elems[38].K2L = deltafamilies[4]
    t0.elems[43].K2L = deltafamilies[5]
    t0.elems[49].K2L = deltafamilies[6]
    t0.elems[74].K2L = deltafamilies[7]
    t0.elems[82].K2L = deltafamilies[8]
    m0 = Map2(t0, terr=None, order=order, nbProc=2)
    #print 'SF1 =',t0.elems[74].K2L, 'and SD0 = ', t0.elems[82].K2L
    sigmaFFS = [sqrt(ex*betx/gamma), sqrt(ex/betx/gamma), sqrt(ey*bety/gamma), sqrt(ey/bety/gamma), 0.01, 0.002]
    h = sqrt(m0.sigma('x', sigmaFFS, gaussian).real)
    j = sqrt(m0.sigma('y', sigmaFFS, gaussian).real)
    h0 = 40e-9 #sigmaFFS[0]
    j0 = 1e-9 #sigmaFFS[2]
    chi = (h-h0)**2/40**2 + (j-j0)**2
    #print 'chi2 = ',chi, 'sig x = ', h, 'sig y = ', j
    return((h-h0)**2/40**2 + (j-j0)**2)
    #return (h)
    #return (j)        
#########################
def banana(deltafamilies):
#########################
    global a, b
    #print "calling banana"
    a = deltafamilies[0]
    b = deltafamilies[1]

    #z = (1-b)**2 + (1-a)**2 
    z = 100*(b-a**2)**2 + (1-a)**2
    #print "try ", deltafamilies, " z =", z
    return (z)

##############
def main():
##############
    global deltaall
    
    deltafamilies1= [-1, 3]
    incr= [4, -4]
    s = Simplex(banana, deltafamilies1, incr)
    values, err, iter = s.minimize()
    print 'args = ', values
    print 'z = ', err
    print 'iterations = ', iter



if __name__ == '__main__':
    main()





