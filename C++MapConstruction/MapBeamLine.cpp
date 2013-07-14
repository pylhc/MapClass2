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

static string drift = "DRIFT";
static string quadrupole = "QUADRUPOLE";
static string sbend = "SBEND";
static string multipole = "MULTIPOLE";
static string decapole = "DECAPOLE";
static string sextupole = "SEXTUPOLE";
static string octupole = "OCTUPOLE";
static string DX = "DX";
static string DY = "DY";


MapBeamLine::MapBeamLine(Twiss t, int order, int nbthreads, int fmultipole) {
	omp_set_num_threads(nbthreads); 
	vector<vector<Polynom<double>>> v = separateComplexList(EQ(4, order));
	Polynom<double> x = X<double>(order);
	Polynom<double> px = PX<double>(order);
	Polynom<double> y = Y<double>(order);
	Polynom<double> py = PY<double>(order);
	Polynom<double> d = D<double>(order);
	Polynom<double> s = S<double>(order);
	Polmap<double> R = generateDefaultMap(x, px, y, py, d, s);
	
	Polmap<double>* Res = new Polmap<double>[nbthreads];
	for (int i = 0; i < nbthreads; i ++)
		Res[i] = generateDefaultMap(x, px, y, py, d, s); 
	int size = t.elems.size();	
	#pragma omp parallel for shared(Res) schedule(static)
	for (int i = 0; i < size; i ++) {
		int index = omp_get_thread_num();
		Polmap<double> mp = mapForElement(t.elems[i], v, x, px, y, py, d, s, fmultipole);	
		Res[index] = mp * Res[index];		
	}
	for (int i = 0; i < nbthreads; i ++)
		R = Res[i] * R;
	polmap = R.getMap();
	for (unordered_map<string, Polynom<double>>:: iterator it = R.pols.begin(); it != R.pols.end(); it ++) 
		pols[it->first] = it->second; 
	delete [] Res;
}

MapBeamLine::MapBeamLine(Twiss t, Twiss terr, int order, int nbthreads, int fmultipole) {
	omp_set_num_threads(nbthreads);
	vector<vector<Polynom<double>>> v = separateComplexList(EQ(4, order));
	Polynom<double> x = X<double>(order);
	Polynom<double> px = PX<double>(order);
	Polynom<double> y = Y<double>(order);
	Polynom<double> py = PY<double>(order);
	Polynom<double> d = D<double>(order);
	Polynom<double> s = S<double>(order);
	Polmap<double> R = generateDefaultMap(x, px, y, py, d, s);
	
	Polmap<double>* Res = new Polmap<double>[nbthreads];
	for (int i = 0; i < nbthreads; i ++)
		Res[i] = generateDefaultMap(x, px, y, py, d, s); 
	
	int size = t.elems.size();	
	#pragma omp parallel for shared(Res) schedule(static)
	for (int i = 0; i < size; i ++) {
		int index = omp_get_thread_num();
		Polmap<double> mp = mapForElement(t.elems[i], v, x, px, y, py, d, s, fmultipole);	
		double dx = atof(terr.elems[i][DX].c_str());
          	double dy = atof(terr.elems[i][DY].c_str());
		mp = mp.eval("x", Polynom<double>(order, 1E-18, "x", 1) + dx); 
		mp = mp.eval("y", Polynom<double>(order, 1E-18, "y", 1) + dy);
		Res[index] = mp * Res[index];		
	}
	for (int i = 0; i < nbthreads; i ++)
		R = Res[i] * R;
	polmap = R.getMap();
	for (unordered_map<string, Polynom<double>>:: iterator it = R.pols.begin(); it != R.pols.end(); it ++) 
		pols[it->first] = it->second; 
	delete [] Res;
}

MapBeamLine::MapBeamLine(string filename, int order, int nbthreads, int fmultipole) {
	Twiss t = Twiss(filename);
	t.stripLine();
	omp_set_num_threads(nbthreads);
	vector<vector<Polynom<double>>> v = separateComplexList(EQ(4, order));
	Polynom<double> x = X<double>(order);
	Polynom<double> px = PX<double>(order);
	Polynom<double> y = Y<double>(order);
	Polynom<double> py = PY<double>(order);
	Polynom<double> d = D<double>(order);
	Polynom<double> s = S<double>(order);
	Polmap<double> R = generateDefaultMap( x, px, y, py, d, s);
	
	Polmap<double>* Res = new Polmap<double>[nbthreads];
	for (int i = 0; i < nbthreads; i ++)
		Res[i] = generateDefaultMap(x, px, y, py, d, s); 
	
	int size = t.elems.size();	
	#pragma omp parallel for shared(Res) schedule(static)
	for (int i = 0; i < size; i ++) {
		int index = omp_get_thread_num();
		Polmap<double> mp = mapForElement(t.elems[i], v, x, px, y, py, d, s, fmultipole);	
		Res[index] = mp * Res[index];		
	}
	for (int i = 0; i < nbthreads; i ++)
		R = Res[i] * R;
	polmap = R.getMap();
	for (unordered_map<string, Polynom<double>>:: iterator it = R.pols.begin(); it != R.pols.end(); it ++) 
		pols[it->first] = it->second; 
	delete [] Res; 
}

MapBeamLine::MapBeamLine(string filename, string filenameerr, int order, int nbthreads, int fmultipole) {
	Twiss t = Twiss(filename);
	t.stripLine();
	Twiss terr = Twiss(filenameerr);
	terr.stripLine();
	omp_set_num_threads(nbthreads);
	vector<vector<Polynom<double>>> v = separateComplexList(EQ(4, order));
	Polynom<double> x = X<double>(order);
	Polynom<double> px = PX<double>(order);
	Polynom<double> y = Y<double>(order);
	Polynom<double> py = PY<double>(order);
	Polynom<double> d = D<double>(order);
	Polynom<double> s = S<double>(order);
	Polmap<double> R = generateDefaultMap(x, px, y, py, d, s);
	
	Polmap<double>* Res = new Polmap<double>[nbthreads];
	for (int i = 0; i < nbthreads; i ++)
		Res[i] = generateDefaultMap(x, px, y, py, d, s); 
	
	int size = t.elems.size();	
	#pragma omp parallel for shared(Res) schedule(static)
	for (int i = 0; i < size; i ++) {
		int index = omp_get_thread_num();
		Polmap<double> mp = mapForElement(t.elems[i], v, x, px, y, py, d, s, fmultipole);	
		double dx = atof(terr.elems[i][DX].c_str());
          	double dy = atof(terr.elems[i][DY].c_str());
		mp = mp.eval("x", Polynom<double>(order, 1E-18, "x", 1) + dx); 
		mp = mp.eval("y", Polynom<double>(order, 1E-18, "y", 1) + dy);
		Res[index] = mp * Res[index];		
	}
	for (int i = 0; i < nbthreads; i ++)
		R = Res[i] * R;
	polmap = R.getMap();
	for (unordered_map<string, Polynom<double>>:: iterator it = R.pols.begin(); it != R.pols.end(); it ++) 
		pols[it->first] = it->second; 
	delete [] Res; 
}



