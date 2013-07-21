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
#include "Twiss.h"

Twiss::Twiss(std::string filename) {
	string line;
  	ifstream file (filename);
  	vector<string> labels;
    	vector<string> types;
    	unordered_map<string, string> e;
  	if (file.is_open())
  	{
    		while (!file.eof())
    		{
    			
      		getline (file,line);
      		if (file.eof())  break;	
      		istringstream iss (line);
      		vector<string> splt;
      		do
    			{
        			string sub;
        			iss >> sub; 
        			splt.push_back(sub);
    			} while (iss);

      		if (line.find("@ ") != string::npos && line.find("%") != string::npos && line.find("s") == string::npos) {
      			types_parameters[splt[1]] = "%le";
      			parameters[splt[1]] = splt[3]; 
      		}
      		else if (line.find("@ ") != string::npos && line.find("%") != string::npos && line.find("s") != string::npos) {
      			trim_if(splt[1],is_any_of(":"));		
      			trim_if(splt[3],is_any_of("\""));
      			string s = splt[3];
      			if (splt.size() > 3) {
      				for (unsigned int i = 4; i < splt.size(); i ++)
      					s += " " + splt[i];
      				trim_if(s,is_any_of("\""));
      			}
      			types_parameters[splt[1]] = "%s";
      			parameters[splt[1]] = s;
      		}
      		if (line.find("* ") != string::npos || line.find("*\t") != string::npos) {
      			labels = splt;
      			labels.erase(labels.begin());
      		}
      		if (line.find("$ ") != string::npos || line.find("$\t") != string::npos) {
      			types = splt;
      			types.erase(types.begin());
      			for (unsigned int i = 0; i < labels.size(); i ++)
      				types_parameters[labels[i]] = types[i];
      		}
      		if (line.find("@") == string::npos && line.find("%") == string::npos && line.find("*") == string::npos) {
        			
        			for (unsigned int j = 0; j < labels.size(); j ++)  {
          		/*	if (types[j].compare("%d") == 0)
          				e[labels[j]] = atoi(splt[j].c_str());
          			else if (types[j].compare("%le") == 0) {
          					e[labels[j]] = atof(splt[j].c_str());
          				}
            				else if (types[j].compare("%s") == 0) {
          						trim_if(splt[j], is_string_of(" "));
          						trim_if(splt[j], is_string_of("\""));
          						e[labels[j]] = splt[j];
          					} 
          			*/

          			trim_if(splt[j], is_any_of("\" "));
          			e[labels[j]] = splt[j];
          			
          		}
          		if (line.find("$") != string::npos) 
          			markers.push_back(e);
          		else {
				
				//string keyword = e["KEYWORD"];
                                //strip the line
				//if (keyword.compare("MARKER") != 0 && keyword.compare("MONITOR") != 0 
			//		&& keyword.compare("MATRIX") != 0)

					elems.push_back(e);
			}	
      		}
      		
    		}
    		file.close();
  	}
  	else cout << "Unable to open file"; 
}

Twiss::Twiss(unordered_map<string, string> p, vector<unordered_map<string, string>> e, vector<unordered_map<string, string>> m, unordered_map<string, string> types){
	parameters = move(p);
	elems = move(e);
	markers = move(m);
	types_parameters = move(types);
}

void Twiss:: printtwiss() {
	for (unordered_map<string, string>::iterator i = parameters.begin(); i != parameters.end(); ++i) {
		cout <<i->first << " " << i->second << endl;
	}
	cout << markers.size();
	for (vector<unordered_map<string, string>>:: iterator i = elems.begin(); i != elems.end(); ++i) {
		unordered_map<string, string> e = *i;
		for (unordered_map<string, string>::iterator j = e.begin(); j != e.end(); ++j) {
			cout <<j->first << " "<<j->second << endl;
		}
	
	}
}

void Twiss:: stripLine() {
	vector<unordered_map<string, string>>:: iterator i = elems.begin();
	while (i != elems.end()) {
		unordered_map<string, string> e = *i;
		string keyword = e["KEYWORD"];
		if (keyword.compare("MARKER") == 0 || keyword.compare("MONITOR") == 0 || keyword.compare("MATRIX") == 0)
                	i = elems.erase(i);
		else ++i;
	}
}

