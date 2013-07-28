#include <iostream>
#include <fstream>
#include <algorithm>
#include <string>
#include <list>
#include <unordered_map>
#include <time.h>
#include "Transport.h"
#include "Twiss.h"

using namespace std;
using namespace boost;

class MapBeamLine : public Polmap<double> {
	public:
		string polmap;
	public:
		MapBeamLine(string filename, int order, int nbthreads, int fmultipole, bool strpl);
		MapBeamLine(string filename, string filenamerr, int order, int nbthreads, int fmultipole, bool strpl);
		MapBeamLine(Twiss t, int order, int nbthreads, int fmultipole, bool strpl);
		MapBeamLine(Twiss t, Twiss terr, int order, int nbthreads, int fmultipole, bool strpl);
};




