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

BOOST_PYTHON_MODULE(mapbeamline_wrapper)
{
     def("constructMapFromTwissFile", _TwissFile_wrapper);
     def("constructMapFromTwissFileWithErr", _TwissFileErr_wrapper);
}
