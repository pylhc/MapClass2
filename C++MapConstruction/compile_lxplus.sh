#!/bin/bash 


export CPLUS_INCLUDE_PATH=../libs/boost_1_53_0/include/
g++ -O2 -fPIC -c Twiss.cpp --std=c++0x  -fopenmp 
g++ -O2 -fPIC -c MapBeamLine.cpp --std=c++0x  -fopenmp
g++ -O2 -fPIC -c MapBeamLine_wrapper.cpp --std=c++0x -fopenmp -I/usr/include/python2.6
g++ -shared -o mapbeamline.so Twiss.o MapBeamLine.o MapBeamLine_wrapper.o ../libs/boost_1_53_0/libboost_python.so.1.53.0 -lgomp
cp mapbeamline.so ../doc/FFSexample
