from ctypes import *
lib = cdll.LoadLibrary('./libMapBeamLine.so')


def MapBeamLine(filename, order, nbthreads):
	return lib.MapBeamLine_new(filename, order, nbthreads)

if __name__ == "__main__":
	lib.MapBeamLine_new.argtypes = [c_char_p, c_int, c_int]
	lib.MapBeamLine_new.restype = c_char_p
	s = lib.MapBeamLine_new("/home/diana/Thesis/cernthesis/MapClass2/doc/FFSexample/assets/ffs.twiss", 3, 2).split("|")
	print s
	
