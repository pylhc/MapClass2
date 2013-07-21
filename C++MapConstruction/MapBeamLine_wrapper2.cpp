#include <vector>
#include "MapBeamLine.h"
#include <Python.h>
#include "boost/python.hpp"
#include "boost/python/module.hpp"
#include "boost/python/def.hpp"
using namespace std;
using namespace boost::python;

str _TwissFile_wrapper(str filename, int order, int nbthreads, int fmultipole, bool strpl) {
	string file = extract<std::string>(filename);
	MapBeamLine mp = MapBeamLine(file, order, nbthreads, fmultipole, strpl);
	return mp.polmap.c_str();
}

str _TwissFileErr_wrapper(str filename, str filenameerr, int order, int nbthreads, int fmultipole, bool strpl) {
	string file = extract<std::string>(filename);
	string fileerr = extract<std::string>(filenameerr);
	MapBeamLine mp = MapBeamLine(file, fileerr, order, nbthreads, fmultipole, strpl);
	return mp.polmap.c_str();
}

str _TwissObject_wrapper(dict tobject, int order, int nbthreads, int fmultipole, bool strpl) {

        boost::python::list elements = extract<boost::python::list>(tobject["elems"]);
        boost::python::list markers = extract<boost::python::list>(tobject["markers"]);
        boost::python::dict types = extract<boost::python::dict>(tobject["types_parameters"]);

        boost::python::list keys = types.keys();
        int len_keys = extract<int>(keys.attr("__len__")());
        unordered_map<string, string> tps;
        for (int j = 0; j < len_keys; ++j) {
                std::string key = extract<std::string> (keys[j]);
                tps[key] = extract<std::string>(types[key]);
        }
        dict py_dict;
        keys = tobject.keys();
        len_keys = extract<int>(keys.attr("__len__")());
        unordered_map<string, string> parameters;
        for (int j = 0; j < len_keys; ++j) {
                std::string key = extract<std::string> (keys[j]);
                if (key.compare("elems") != 0 && key.compare("markers") != 0 && key.compare("types_parameters") != 0) {
			if (tps[key].compare("%s") == 0)
                       		parameters[key] = extract<std::string>(tobject[key]);
                	else
                        	parameters[key] = boost::lexical_cast<std::string>(extract<double>(tobject[key]));
		}

        }

        int elements_length = extract<int>(elements.attr("__len__")());
        vector<unordered_map<string, string>> elems;
        unordered_map<string, string> map;
        for (int i = 0; i < elements_length; i++) {
                py_dict = extract<dict>(elements[i]);
                keys = py_dict.keys();
                len_keys = extract<int>(keys.attr("__len__")());
                for (int j = 0; j < len_keys; ++j) {
                        std::string key = boost::python::extract<std::string> (keys[j]);
                       if (tps[key].compare("%s") == 0)
                                map[key] = extract<std::string>(py_dict[key]);
                        else
                                map[key] = boost::lexical_cast<std::string>(extract<double>(py_dict[key]));
                }
                elems.push_back(map);
        }

        int markers_length = extract<int>(markers.attr("__len__")());
        vector<unordered_map<string, string>> mks;
        for (int i = 0; i < markers_length; i++) {
                py_dict = extract<dict>(markers[i]);
                keys = py_dict.keys();
                len_keys = extract<int>(keys.attr("__len__")());
                for (int j = 0; j < len_keys; ++j) {
                        std::string key = boost::python::extract<std::string> (keys[j]);
                        if (tps[key].compare("%s") == 0)
                                map[key] = extract<std::string>(py_dict[key]);
                        else
                                map[key] = boost::lexical_cast<std::string>(extract<double>(py_dict[key]));
                }
                mks.push_back(map);
        }

        Twiss t = Twiss(parameters, elems, mks, tps);
        MapBeamLine mp = MapBeamLine(t, order, nbthreads, fmultipole, strpl);
        return mp.polmap.c_str();
}

str _TwissObjectErr_wrapper(dict tobject, dict terrobject, int order, int nbthreads, int fmultipole, bool strpl) {
        boost::python::list elements = extract<boost::python::list>(tobject["elems"]);
        boost::python::list markers = extract<boost::python::list>(tobject["markers"]);
        boost::python::dict types = extract<boost::python::dict>(tobject["types_parameters"]);

        boost::python::list keys = types.keys();
        int len_keys = extract<int>(keys.attr("__len__")());
        unordered_map<string, string> tps;
       	for (int j = 0; j < len_keys; ++j) {
                std::string key = extract<std::string> (keys[j]);
                tps[key] = extract<std::string>(types[key]);
        }
        dict py_dict;
        keys = tobject.keys();
        len_keys = extract<int>(keys.attr("__len__")());
        unordered_map<string, string> parameters;
        for (int j = 0; j < len_keys; ++j) {
                std::string key = extract<std::string> (keys[j]);
                 if (key.compare("elems") != 0 && key.compare("markers") != 0 && key.compare("types_parameters") != 0) {
			if (tps[key].compare("%s") == 0)
                        	parameters[key] = extract<std::string>(tobject[key]);
                	else
                        	parameters[key] = boost::lexical_cast<std::string>(extract<double>(tobject[key]));
		}

        }

        int elements_length = extract<int>(elements.attr("__len__")());
        vector<unordered_map<string, string>> elems;
        unordered_map<string, string> map;
        for (int i = 0; i < elements_length; i++) {
                py_dict = extract<dict>(elements[i]);
                keys = py_dict.keys();
                len_keys = extract<int>(keys.attr("__len__")());
                for (int j = 0; j < len_keys; ++j) {
                        std::string key = boost::python::extract<std::string> (keys[j]);
                        if (tps[key].compare("%s") == 0)
                                map[key] = extract<std::string>(py_dict[key]);
                        else
                                map[key] = boost::lexical_cast<std::string>(extract<double>(py_dict[key]));
                }
                elems.push_back(map);
        }

        int markers_length = extract<int>(markers.attr("__len__")());
        vector<unordered_map<string, string>> mks;
        for (int i = 0; i < markers_length; i++) {
                py_dict = extract<dict>(markers[i]);
                keys = py_dict.keys();
                len_keys = extract<int>(keys.attr("__len__")());
                for (int j = 0; j < len_keys; ++j) {
                        std::string key = boost::python::extract<std::string> (keys[j]);
                        if (tps[key].compare("%s") == 0)
                                map[key] = extract<std::string>(py_dict[key]);
                        else
                                map[key] = boost::lexical_cast<std::string>(extract<double>(py_dict[key]));
        	}
                mks.push_back(map);
        }

        Twiss t = Twiss(parameters, elems, mks, tps);

        elements = extract<boost::python::list>(terrobject["elems"]);
        markers = extract<boost::python::list>(terrobject["markers"]);
        types = extract<boost::python::dict>(terrobject["types_parameters"]);

        terrobject["elems"].del();
        terrobject["markers"].del();
        terrobject["types_parameters"].del();

        keys = types.keys();
        len_keys = extract<int>(keys.attr("__len__")());
        unordered_map<string, string> tpserr;
        for (int j = 0; j < len_keys; ++j) {
                std::string key = extract<std::string> (keys[j]);
                tpserr[key] = extract<std::string>(types[key]);
        }
        keys = terrobject.keys();
        len_keys = extract<int>(keys.attr("__len__")());
        unordered_map<string, string> parameterserr;
        for (int j = 0; j < len_keys; ++j) {
                std::string key = extract<std::string> (keys[j]);
                if (tps[key].compare("%s") == 0)
                        parameterserr[key] = extract<std::string>(terrobject[key]);
                else
                        parameterserr[key] = boost::lexical_cast<std::string>(extract<double>(terrobject[key]));

        }

        elements_length = extract<int>(elements.attr("__len__")());
        vector<unordered_map<string, string>> elemserr(elements_length);
        for (int i = 0; i < elements_length; i++) {
                py_dict = extract<dict>(elements[i]);
               keys = py_dict.keys();
                len_keys = extract<int>(keys.attr("__len__")());
                for (int j = 0; j < len_keys; ++j) {
                        std::string key = boost::python::extract<std::string> (keys[j]);
                        if (tpserr[key].compare("%s") == 0)
                                map[key] = extract<std::string>(py_dict[key]);
                        else
                                map[key] = boost::lexical_cast<std::string>(extract<double>(py_dict[key]));
                }
                elemserr.push_back(map);
        }

        markers_length = extract<int>(markers.attr("__len__")());
        vector<unordered_map<string, string>> mkserr(markers_length);
        for (int i = 0; i < markers_length; i++) {
                py_dict = extract<dict>(markers[i]);
                keys = py_dict.keys();
                len_keys = extract<int>(keys.attr("__len__")());
                for (int j = 0; j < len_keys; ++j) {
                        std::string key = boost::python::extract<std::string> (keys[j]);
                        if (tpserr[key].compare("%s") == 0)
                                map[key] = extract<std::string>(py_dict[key]);
                        else
                                map[key] = boost::lexical_cast<std::string>(extract<double>(py_dict[key]));
                }
                mkserr.push_back(map);
        }

        Twiss terr = Twiss(parameterserr, elemserr, mkserr, tpserr);
        MapBeamLine mp = MapBeamLine(t, terr, order, nbthreads, fmultipole, strpl);
        return mp.polmap.c_str();
}

BOOST_PYTHON_MODULE(mapbeamline)
{
     def("constructMapFromTwissFile2", _TwissFile_wrapper);
     def("constructMapFromTwissFileWithErr2", _TwissFileErr_wrapper);
     def("constructMapFromTwissObject2", _TwissObject_wrapper);
     def("constructMapFromTwissObjectWithErr2", _TwissObjectErr_wrapper);
}
