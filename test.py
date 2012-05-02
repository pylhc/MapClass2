#!/usr/bin/env python

import math
import argparse

from mapclass import *

parser = argparse.ArgumentParser(description='Tester for MapClass')
parser.add_argument('-o', help='Run the test at a given order',
                    action='store', type=int, dest='order', default=6)
args = parser.parse_args()

betx=66.14532014
bety=17.92472388
gamma=3e6
ex=68e-8
ey=2e-8
sigmaFFS=[sqrt(ex*betx/gamma), sqrt(ex/betx/gamma), sqrt(ey*bety/gamma), sqrt(ey/bety/gamma), 0.01]

file='assets/fort.18'
map=Map2(args.order,file)
print "sigmax=",sqrt(map.sigma('fx',sigmaFFS)),";"
print "sigmay=",sqrt(map.sigma('fy',sigmaFFS)),";"
print "sigmapx=",sqrt(map.sigma('fpx',sigmaFFS)),";"
print "sigmapy=",sqrt(map.sigma('fpy',sigmaFFS)),";"
print "offsetx=",map.offset('fx',sigmaFFS),";"
print "offsety=",map.offset('fy',sigmaFFS),";"
print "offsetpx=",map.offset('fpx',sigmaFFS),";"
print "offsetpy=",map.offset('fpy',sigmaFFS),";"

print "\n\n########## GaussianDelta ##########"

map=Map2(args.order,file)
print "sigmax=",sqrt(map.sigma('fx',sigmaFFS, True)),";"
print "sigmay=",sqrt(map.sigma('fy',sigmaFFS, True)),";"
print "sigmapx=",sqrt(map.sigma('fpx',sigmaFFS, True)),";"
print "sigmapy=",sqrt(map.sigma('fpy',sigmaFFS, True)),";"
print "offsetx=",map.offset('fx',sigmaFFS, True),";"
print "offsety=",map.offset('fy',sigmaFFS, True),";"
print "offsetpx=",map.offset('fpx',sigmaFFS, True),";"
print "offsetpy=",map.offset('fpy',sigmaFFS, True),";"
