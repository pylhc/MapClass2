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
#include "Twiss.h"

Twiss::Twiss(std::string filename) {
	string line;
  	ifstream file (filename);
  	vector<string> labels;
    vector<string> types;
    unordered_map<string, any> e;
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
      			parameters[splt[1]] = atof(splt[3].c_str()); 
      		}
      		else if (line.find("@ ") != string::npos && line.find("%") != string::npos && line.find("s") != string::npos) {
      			trim_if(splt[1],is_any_of(":"));		
      			trim_if(splt[3],is_any_of("\""));
      			string s = splt[3];
      			if (splt.size() > 3) {
      				for (int i = 4; i < splt.size(); i ++)
      					s += " " + splt[i];
      				trim_if(s,is_any_of("\""));
      			}
      			parameters[splt[1]] = s;
      		}
      		if (line.find("* ") != string::npos || line.find("*\t") != string::npos) {
      			labels = splt;
      			labels.erase(labels.begin());
      		}
      		if (line.find("$ ") != string::npos || line.find("$\t") != string::npos) {
      			types = splt;
      			types.erase(types.begin());
      		}
      		if (line.find("@") == string::npos && line.find("%") == string::npos && line.find("*") == string::npos) {
        			e.clear();
        			
        			for (int j = 0; j < labels.size(); j ++)  {
          			if (types[j].compare("%hd") == 0)
          				e[labels[j]] = atoi(splt[j].c_str());
          			else if (types[j].compare("%le") == 0) {
          					e[labels[j]] = atof(splt[j].c_str());
          				}
            				else if (types[j].compare("%s") == 0) {
          						trim_if(splt[j], is_any_of(" "));
          						trim_if(splt[j], is_any_of("\""));
          						e[labels[j]] = splt[j];
          					} 
          		}
          		if (line.find("$") != string::npos) 
          			markers.push_back(e);
          		else elems.push_back(e);	
      		}
      		
    		}
    		file.close();
  	}
  	else cout << "Unable to open file"; 
}

Twiss::Twiss(unordered_map<string, any> p, vector<unordered_map<string, any>> e, vector<unordered_map<string, any>> m){
	parameters = move(p);
	elems = move(e);
	markers = move(m);
}

void Twiss:: printtwiss() {
	for (unordered_map<string, any>::iterator i = parameters.begin(); i != parameters.end(); ++i) {
		cout <<i->first << " ";
		try
    		{
        		cout << any_cast<double>(i->second)<<endl;
        	}
    		catch(const boost::bad_any_cast &)
   	 	{
        		cout << any_cast<string>(i->second)<<endl;
    		}	
	}
	cout << markers.size();
	for (vector<unordered_map<string, any>>:: iterator i = elems.begin(); i != elems.end(); ++i) {
		unordered_map<string, any> e = *i;
		for (unordered_map<string, any>::iterator j = e.begin(); j != e.end(); ++j) {
			cout <<j->first << " ";
			try
    			{
        			cout << any_cast<double>(j->second)<<endl;
        		}
    			catch(const boost::bad_any_cast &)
   	 		{
        			try
    				{
        				cout << any_cast<string>(j->second)<<endl;
        			}
        			catch(const boost::bad_any_cast &)
   	 			{
        				cout << any_cast<int>(j->second)<<endl;
    				}
    			}	
		}
	
	}
}

