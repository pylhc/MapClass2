#include <iostream>
#include <fstream>
#include <algorithm>
#include <list>
#include <string>
#include <cstdlib>
#include <vector>
#include <unordered_map>
#include <boost/any.hpp>
#include <boost/algorithm/string.hpp>
#include <boost/algorithm/string/trim.hpp>
#include <sstream>

#ifndef TWISS_H
#define TWISS_H

using namespace std;
using namespace boost;

class Twiss {
	public:
		unordered_map<string, any> parameters;
		vector<unordered_map<string, any>> elems;
		vector<unordered_map<string, any>> markers;
	public:
		Twiss(std::string s);
		Twiss(unordered_map<string, any> p, vector<unordered_map<string, any>> e, vector<unordered_map<string, any>> m);
		void printtwiss();
};


#endif
