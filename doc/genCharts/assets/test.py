#! /usr/bin/env python

import sys

if len(sys.argv) < 2:
  sys.stderr.write('Usage: sys.argv[0] dir')
  sys.exit(1)

sys.path.append(sys.argv[1])

from metaclass import twiss
from mapclass import Map2

t = twiss()

m = Map2(t)

mm = Map2(order=6)

v = ['x','px','y','py']

print mm.compc(m,v).real
