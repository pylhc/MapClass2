#include <iostream>
#include <fstream>
#include <algorithm>
#include <list>
#include <string>
#include <cstdlib>
#include <vector>
#include <unordered_map>
#include "boost/algorithm/string.hpp"
#include "boost/algorithm/string/trim.hpp"
#include <sstream>

#ifndef TWISS_H
#define TWISS_H

using namespace std;
using namespace boost;

class Twiss {
	public:
		unordered_map<string, string> parameters;
		vector<unordered_map<string, string>> elems;
		vector<unordered_map<string, string>> markers;
		unordered_map<string, string> types_parameters;
	public:
		Twiss(std::string s);
		Twiss(unordered_map<string, string> p, vector<unordered_map<string, string>> e, vector<unordered_map<string, string>> m, unordered_map<string, string> t);
		void printtwiss();
		void stripLine();
};


#endif
