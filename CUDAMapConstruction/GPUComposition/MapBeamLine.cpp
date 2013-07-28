#include <iostream>
#include <fstream>
#include <algorithm>
#include <string>
#include <list>
#include <unordered_map>
#include "MapBeamLine.h"

#define SIZE	4096

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

MapBeamLine::MapBeamLine(Twiss t, int order, int nbthreads) {
	cudaError_t error;
	int devID = 0;
    	error = cudaGetDevice(&devID);

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
		Polmap<double> mp = mapForElement(t.elems[i], v, x, px, y, py, d, s);	
		//if (mp.pols.size() != 0)
			Res[index] = mp * Res[index];		
	}

	R = Res[0];
	double start = omp_get_wtime();
	for (int i = 1; i < nbthreads; i ++)
	//   	R = Res[i] * R;
		R = composeGPU(Res[i], R);
	//	R = Res[i].compositionWithoutHash(R);
	double end = omp_get_wtime();
	cout << 1000 * (end - start) << " ms" << endl; 

//	polmap = R.getMap();
	for (unordered_map<string, Polynom<double>>:: iterator it = R.pols.begin(); it != R.pols.end(); it ++) 
	{	pols[it->first] = it->second;
// 		cout << it->second.terms.size() << endl;
	}
	delete [] Res;
	cudaDeviceReset();
}

MapBeamLine::MapBeamLine(Twiss t, Twiss terr, int order, int nbthreads) {
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
		Polmap<double> mp = mapForElement(t.elems[i], v, x, px, y, py, d, s);	
		double dx = atof(terr.elems[i][DX].c_str());
          	double dy = atof(terr.elems[i][DY].c_str());
		mp = mp.eval("x", Polynom<double>(order, 1E-18, "x", 1) + dx); 
		mp = mp.eval("y", Polynom<double>(order, 1E-18, "y", 1) + dy);
		//if (mp.pols.size() != 0)
			Res[index] = mp * Res[index];		
	}
	for (int i = 0; i < nbthreads; i ++)
		R = Res[i] * R;
	polmap = R.getMap();
	//for (unordered_map<string, Polynom<double>>:: iterator it = R.pols.begin(); it != R.pols.end(); it ++) 
	//	pols[it->first] = it->second; 
	delete [] Res;
}

MapBeamLine::MapBeamLine(string filename, int order, int nbthreads) {
	Twiss t = Twiss(filename);
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
		Polmap<double> mp = mapForElement(t.elems[i], v, x, px, y, py, d, s);	
		//if (mp.pols.size() != 0)
			Res[index] = mp * Res[index];		
	}
	for (int i = 0; i < nbthreads; i ++)
		R = Res[i] * R;
	polmap = R.getMap();
	//for (unordered_map<string, Polynom<double>>:: iterator it = R.pols.begin(); it != R.pols.end(); it ++) 
	//	pols[it->first] = it->second; 
	delete [] Res; 
}

MapBeamLine::MapBeamLine(string filename, string filenameerr, int order, int nbthreads) {
	Twiss t = Twiss(filename);
	Twiss terr = Twiss(filenameerr);
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
		Polmap<double> mp = mapForElement(t.elems[i], v, x, px, y, py, d, s);	
		double dx = atof(terr.elems[i][DX].c_str());
          	double dy = atof(terr.elems[i][DY].c_str());
		mp = mp.eval("x", Polynom<double>(order, 1E-18, "x", 1) + dx); 
		mp = mp.eval("y", Polynom<double>(order, 1E-18, "y", 1) + dy);
		//if (mp.pols.size() != 0)
			Res[index] = mp * Res[index];		
	}
	for (int i = 0; i < nbthreads; i ++)
		R = Res[i] * R;
	polmap = R.getMap();
	//for (unordered_map<string, Polynom<double>>:: iterator it = R.pols.begin(); it != R.pols.end(); it ++) 
	//	pols[it->first] = it->second; 
	delete [] Res; 
}

int main(int argc,char *argv[]) {
        int order = 3;
        int nbthreads = 1;
        string filename = "";
        for (int i = 1; i < argc; i += 2) {
                if (strcmp(argv[i], "-o") == 0)
                        order = atoi(argv[i + 1]);
                else if (strcmp(argv[i], "-n") == 0)
                                nbthreads = atoi(argv[i + 1]);
                else if (strcmp(argv[i], "-f") == 0)
                                filename = argv[i + 1];
        }
        Twiss t = Twiss(filename);

        double start = omp_get_wtime();
        MapBeamLine mbl = MapBeamLine(t, order, nbthreads);
        double end = omp_get_wtime();
        cout << "running time = " <<  end - start << endl;
        //mbl.printpolmap();
        return 0;
}


