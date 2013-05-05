#!/bin/bash 

sudo apt-get install libboost-all-devs

g++ -O2 -fPIC -c Twiss.cpp --std=c++0x
g++ -O2 -fPIC -c MapBeamLine.cpp --std=c++0x -fopenmp 
g++ -O2 -fPIC -c MapBeamLine_wrapper.cpp --std=c++0x  -lboost_python -I /usr/include/python2.7 -fopenmp
g++ -shared -o mapbeamline_wrapper.so Twiss.o MapBeamLine.o MapBeamLine_wrapper.o -lboost_python -lgomp
cp mapbeamline_wrapper.so ../doc/FFSexample
