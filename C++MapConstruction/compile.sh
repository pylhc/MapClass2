#!/bin/bash 

source /afs/cern.ch/sw/lcg/contrib/gcc/4.4.3/x86_64-slc5/setup.sh
export CPLUS_INCLUDE_PATH=../libs/boost_1_53_0/include/
g++ -O2 -fPIC -c Twiss.cpp --std=c++0x  -fopenmp 
g++ -O2 -fPIC -c MapBeamLine.cpp --std=c++0x  -fopenmp
g++ -O2 -fPIC -c MapBeamLine_wrapper.cpp --std=c++0x -fopenmp -I/usr/include/python2.6
g++ -shared -o mapbeamline_wrapper.so Twiss.o MapBeamLine.o MapBeamLine_wrapper.o ../libs/boost_1_53_0/lib/libboost_python.so.1.53.0 -lgomp
cp mapbeamline_wrapper.so ../doc/FFSexample
