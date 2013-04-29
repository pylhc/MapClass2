#include <iostream>
#include <fstream>
#include <algorithm>
#include <string>
#include <list>
#include <unordered_map>
#include <time.h>
#include "MapBeamLine.h"

using namespace std;
using namespace boost;

timespec diff(timespec start, timespec end);
static string drift = "DRIFT";
static string quadrupole = "QUADRUPOLE";
static string sbend = "SBEND";
static string multipole = "MULTIPOLE";
static string decapole = "DECAPOLE";
static string sextupole = "SEXTUPOLE";
static string octupole = "OCTUPOLE";

void MapBeamLine::fromTwiss(Twiss t, Twiss terr, int order, int nbthreads) {
	omp_set_num_threads(nbthreads);
	timespec time1, time2; 
	vector<vector<Polynom<double>>> v = separateComplexList(EQ(4, order));
	Polynom<double> x = X<double>(order);
	Polynom<double> px = PX<double>(order);
	Polynom<double> y = Y<double>(order);
	Polynom<double> py = PY<double>(order);
	Polynom<double> d = D<double>(order);
	Polynom<double> s = S<double>(order);
	Polmap<double> R = generateDefaultMap(order, x, px, y, py, d, s);
	
	Polmap<double>* Res = new Polmap<double>[nbthreads];
	for (int i = 0; i < nbthreads; i ++)
		Res[i] = generateDefaultMap(order, x, px, y, py, d, s); 
	
	int size = t.elems.size();	
	#pragma omp parallel for shared(Res) schedule(static)
	for (int i = 0; i < size; i ++) {
		int index = omp_get_thread_num();
		Polmap<double> mp = mapForElement(t.elems[i], order, v, x, px, y, py, d, s);	
		if (mp.pols.size() != 0)
				Res[index] = mp * Res[index];		
	}
	for (int i = 0; i < nbthreads; i ++)
		R = Res[i] * R;
	for (unordered_map<string, Polynom<double>>:: iterator it = R.pols.begin(); it != R.pols.end(); it ++) 
		pols[it->first] = it->second; 
	delete [] Res;
}

MapBeamLine::MapBeamLine(string filename, int order, int nbthreads) {
	Twiss t = Twiss(filename);
	omp_set_num_threads(nbthreads);
	timespec time1, time2; 
	vector<vector<Polynom<double>>> v = separateComplexList(EQ(4, order));
	Polynom<double> x = X<double>(order);
	Polynom<double> px = PX<double>(order);
	Polynom<double> y = Y<double>(order);
	Polynom<double> py = PY<double>(order);
	Polynom<double> d = D<double>(order);
	Polynom<double> s = S<double>(order);
	Polmap<double> R = generateDefaultMap(order, x, px, y, py, d, s);
	
	Polmap<double>* Res = new Polmap<double>[nbthreads];
	for (int i = 0; i < nbthreads; i ++)
		Res[i] = generateDefaultMap(order, x, px, y, py, d, s); 
	
	int size = t.elems.size();	
	#pragma omp parallel for shared(Res) schedule(static)
	for (int i = 0; i < size; i ++) {
		int index = omp_get_thread_num();
		Polmap<double> mp = mapForElement(t.elems[i], order, v, x, px, y, py, d, s);	
		if (mp.pols.size() != 0)
				Res[index] = mp * Res[index];		
	}
	for (int i = 0; i < nbthreads; i ++)
		R = Res[i] * R;
	polmap = R.getMap();
	for (unordered_map<string, Polynom<double>>:: iterator it = R.pols.begin(); it != R.pols.end(); it ++) 
		pols[it->first] = it->second; 
	delete [] Res; 
}

timespec diff(timespec start, timespec end)
{
	timespec temp;
	if ((end.tv_nsec-start.tv_nsec)<0) {
		temp.tv_sec = end.tv_sec-start.tv_sec-1;
		temp.tv_nsec = 1000000000+end.tv_nsec-start.tv_nsec;
	} else {
		temp.tv_sec = end.tv_sec-start.tv_sec;
		temp.tv_nsec = end.tv_nsec-start.tv_nsec;
	}
	return temp;
}


extern "C" {
    const char* MapBeamLine_new(char* filename, int order, int nbthreads){ 
		MapBeamLine mp = MapBeamLine(filename, order, nbthreads);
		return mp.polmap.c_str();
 
	}
    void MapBeamLine_fromTwiss(MapBeamLine* mb, Twiss t, Twiss err, int order, int nbthreads){ 
		mb->fromTwiss(t, err, order, nbthreads); 
	}
}


