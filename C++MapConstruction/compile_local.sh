#!/bin/bash 

g++ -O2 -fPIC -c Twiss.cpp --std=c++0x
g++ -O2 -fPIC -c MapBeamLine.cpp --std=c++0x -fopenmp 
g++ -O2 -fPIC -c MapBeamLine_wrapper.cpp --std=c++0x  -lboost_python -I /usr/include/python2.7 -fopenmp
g++ -shared -o mapbeamline.so Twiss.o MapBeamLine.o MapBeamLine_wrapper.o -lboost_python -lgomp
cp mapbeamline.so ../doc/FFSexample
