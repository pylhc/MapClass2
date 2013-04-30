#include <vector>
#include "MapBeamLine.h"
#include <Python.h>
#include <boost/python.hpp>
#include <boost/python/module.hpp>
#include <boost/python/def.hpp>
using namespace std;
using namespace boost::python;


str _TwissFile_wrapper(str filename, int order, int nbthreads) {
	MapBeamLine mp = MapBeamLine(extract<std::string>(filename), order, nbthreads);
	return mp.polmap.c_str();
}


str _TwissFileErr_wrapper(str filename, str filenameerr, int order, int nbthreads) {
	MapBeamLine mp = MapBeamLine(extract<std::string>(filename), extract<std::string>(filenameerr), order, nbthreads);
	return mp.polmap.c_str();
}

str _TwissObject_wrapper(dict tobject, boost::python::list elements, boost::python::list markers, int order, int nbthreads) {

	boost::python::list keys = tobject.keys(); 
	dict py_dict;
	unordered_map<string, any> parameters;
	
	for (int j = 0; j < len(keys); ++j) {  
		std::string key = boost::python::extract<std::string> (keys[j]); 
		parameters[key] = tobject[key];
     }  
	int elements_length = extract<int>(elements.attr("__len__")());
	vector<unordered_map<string, any>> elems(elements_length);
	unordered_map<string, any> map;
	for (int i = 0; i < elements_length; i++) {
		py_dict = extract<dict>(elements[i]);
		keys = py_dict.keys();  
        	for (int j = 0; j < len(keys); ++j) {  
			std::string key = boost::python::extract<std::string> (keys[j]); 
			map[key] = py_dict[key];
     	}  
      	elems.push_back(map);
	}
	
	int markers_length = extract<int>(markers.attr("__len__")());
	vector<unordered_map<string, any>> mks(markers_length);
	for (int i = 0; i < markers_length; i++) {
		py_dict = extract<dict>(markers[i]);
		keys = py_dict.keys();  
        	for (int j = 0; j < len(keys); ++j) {  
			std::string key = boost::python::extract<std::string> (keys[j]); 
			map[key] = py_dict[key];
     	}  
      	mks.push_back(map);
	}
	
	Twiss t = Twiss(parameters, elems, mks);
	MapBeamLine mp = MapBeamLine(t, order, nbthreads);
	return mp.polmap.c_str();
}

str _TwissObjectErr_wrapper(dict tobject, boost::python::list elements, boost::python::list markers, 
					dict terrobject, boost::python::list errelements, boost::python::list errmarkers,
					int order, int nbthreads) {
	boost::python::list keys = tobject.keys(); 
	dict py_dict;
	unordered_map<string, any> parameters;
	for (int j = 0; j < len(keys); ++j) {  
		std::string key = boost::python::extract<std::string> (keys[j]); 
          parameters[key] = tobject[key];
     }  
	int elements_length = extract<int>(elements.attr("__len__")());
	vector<unordered_map<string, any>> elems(elements_length);
	unordered_map<string, any> map;
	for (int i = 0; i < elements_length; i++) {
		py_dict = extract<dict>(elements[i]);
		keys = py_dict.keys();  
        	for (int j = 0; j < len(keys); ++j) {  
			std::string key = boost::python::extract<std::string> (keys[j]); 
			map[key] = py_dict[key];
      	}  
      	elems.push_back(map);
	}
	
	int markers_length = extract<int>(markers.attr("__len__")());
	vector<unordered_map<string, any>> mks(markers_length);
	for (int i = 0; i < markers_length; i++) {
		py_dict = extract<dict>(markers[i]);
		keys = py_dict.keys();  
        	for (int j = 0; j < len(keys); ++j) {  
			std::string key = boost::python::extract<std::string> (keys[j]); 
			map[key] = py_dict[key];
      	}  
      	mks.push_back(map);
	}
	
	Twiss t = Twiss(parameters, elems, mks);
	
	keys = terrobject.keys(); 
	unordered_map<string, any> errparameters;
	for (int j = 0; j < len(keys); ++j) {  
		std::string key = boost::python::extract<std::string> (keys[j]); 
          map[key] = terrobject[key];
     }  
	elements_length = extract<int>(errelements.attr("__len__")());
	vector<unordered_map<string, any>> errelems(elements_length);
	for (int i = 0; i < elements_length; i++) {
		py_dict = extract<dict>(errelements[i]);
		keys = py_dict.keys();  
        	for (int j = 0; j < len(keys); ++j) {  
			std::string key = boost::python::extract<std::string> (keys[j]);
           	map[key] = py_dict[key];
      	}  
      	errelems.push_back(map);
	}
	
	markers_length = extract<int>(errmarkers.attr("__len__")());
	vector<unordered_map<string, any>> errmks(markers_length);
	for (int i = 0; i < markers_length; i++) {
		py_dict = extract<dict>(errmarkers[i]);
		keys = py_dict.keys();  
        	for (int j = 0; j < len(keys); ++j) {  
			std::string key = boost::python::extract<std::string> (keys[j]);
			map[key] = py_dict[key];
      	}  
      	errmks.push_back(map);
	}
	
	Twiss terr = Twiss(errparameters, errelems, errmks);
	MapBeamLine mp = MapBeamLine(t, terr, order, nbthreads);
	return mp.polmap.c_str();
									
}

BOOST_PYTHON_MODULE(mapbeamline_wrapper)
{
     def("constructMapFromTwissFile", _TwissFile_wrapper);
     def("constructMapFromTwissFileWithErr", _TwissFileErr_wrapper);
     def("constructMapFromTwissObject", _TwissObject_wrapper);
     def("constructMapFromTwissObjectWithErr", _TwissObjectErr_wrapper);
}
