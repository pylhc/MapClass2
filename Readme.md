MAPCLASS2
=========
### A code to aid particle accelerator lattice design and optimization

This Readme.md contains a brief description of the code purpose, installation instructions and usage.

For a more comprehensive description, refer to the bibliography section, where you can find :

[1] a conceptual review,

[2] use and detailed python code explanations,

[3] the implemention with C++ and CUDA libraries,

[4] effects of radiation in dipoles and quadrupoles (Oide effect),

[5] the fort.18 file details, in Appendix B,

[6] the MAD-X project,

[7] and the initial MAPCLASS code.

DESCRIPTION
-----------
The transport of a particle beam through out an accelerator requires the calculation of the beam properties, e.g. beamsize in the transverse planes, correlation, etc. Typically, it is done by tracking several hundreds to thousands of particles, with an statistical post-analysis. Such process could be computationally intensive and too slow for very long or complicated transfer lines. On top of that, during the lattice design and iterative optimization one would like to explore several possibilities, and it would be advantageous to reduce the time of calculation.

MAPCLASS was written to perform a fast calculation of the outgoing beam size, given a transport lattice and a beamsize at the input. Therefore, it avoids the time consuming single particle tracking. It has been extended in MapClass2 to calculate twiss functions inside elements, correlations in the phase space and many others.

The main input to MapClass2 is the lattice, either from
 * a fort.18 file generated in MAD-X PTC, or, alternatively
 * a twiss file from MAD-X
 
Naturally, in order to give fast results we must make a simplification. MAPCLASS assumes the beam at the entry point is gaussian in the horizontal phase space (x, px), vertical phase space (y, py) and longitudinal position (t). The beam energy spread distribution could be chosen to be gaussian or uniform. Also, the lattice twiss parameters alfx and alfy must be zero at the entry point.


INSTALL INSTRUCTIONS
--------------------
Two options are available: C++ libraries (recommended),and CUDA libraries (parallelization)

1. From a bash terminal, move to desired location, and then execute the following lines:
```bash
git clone https://github.com/pylhc/MapClass2.git
cd MapClass2
git submodule init
git submodule update
```

2. The Python variables need to be modified. Here are two ways to do it :

* Include the path inside the python script where you call MapClass2, e.g.
```python
### my python script loading MapClass2_C++
import sys, ctypes
sys.path.append('/home/o/MapClass2/')
sys.path.append('/home/o/MapClass2/C++MapConstruction/')
ctypes.cdll.LoadLibrary("/home/o/MapClass2/libs/boost_1_53_0/libboost_python.so.1.53.0")
```

* by adding to the paths to your ~/.bashrc file, e.g.
```bash
### Mapclass2_C++
export PYTHONPATH="/home/o/MapClass2:$PYTHONPATH"
export PYTHONPATH="$PYTHONPATH:/home/o/MapClass2/C++MapConstruction"
export LD_LIBRARY_PATH="$LD_LIBRARY_PATH:/home/o/MapClass2/libs/boost_1_53_0"

```
#### NOTE: In the previous example we show the usage of C++ LIBRARIES (recommended) in the "C++MapConstruction" dir. If you choose CUDA LIBRARIES (it must be build), use the "CUDAMapConstruction" dir

#### COMMENTS: If you are new to git commands  

```bash
git pull #To update
git reset --hard #To recover in an accidental file remove or modification
```

SOFTWARE REQUIREMENTS
---------------------
+ Python 2.6 or higher is required (not tested with Python 3)
   +  numpy library must be installed or added to the 'libs' folder (tested 1.13.3)
+ gcc 4.4.3 or higher.
   +  If working on the cern afs and your gcc is lower than 4.4.3, then do
```bash
          source /afs/cern.ch/sw/lcg/external/gcc/4.4/x86_64-slc5/setup.sh
```
+ (OPTIONAL) CUDA libraries


HOW TO USE
----------
From any Python terminal environment or script, it is possible to load the modules mapclass and metaclass2.
```python
import mapclass
import metaclass2
```
now, a twiss or fort.18 file (from MAD-X 5.X.X) is required to perform beamsize
calculations.

### INPUT FILES
+ fort.18 files from MAD-X PTC

  Check the PTC module docs in MAD-X as it is very powerful. Here is a common MAD-X coding example to generate fort.18 file.
```
      !!! MAD-X PTC Code to produce fort.18
        ptc_create_universe;
        ptc_create_layout,model=2,method=6,nst=10;
        ptc_normal,icase=5,no=8,deltap=0.00;
        ptc_end; 
      !!! End MAD-X Code
```

+ twiss files from MAD-X Twiss command should be generated with at least the following columns:
```
NAME, KEYWORD, S, L, BETX, BETY, ALFX, ALFY, MUX, MUY, DX, DPX, DY, DPY, ANGLE, K1L, K2L, K3L, K4L
```
  ***WARNING: Some MAPCLASS2 functions assume elements referenced to its exit side.***


EXAMPLES
--------
Some examples of use can be found in ./doc/FFSexample


QUICK Q/A
---------
### Using twiss files
  + How to load a twiss file?

    In a python script or shell:
``` python
    import metaclass2
    tw = metaclass2.twiss2("twissfilename")
```
  + How to check what was loaded?

    `len(tw.elems)` tells you how many elements were loaded, while
    `tw.elems[0]` is the first loaded element.

    `tw.elems[0].KEYWORD` is the keyword used in MAD-X for the element zero. MapClass2 uses the same 
    twiss file column names except
    for `TNAME` which is `TableNAME` (TNAME does not exist in MAD-X).

    `tw.elems[0].BETX` is the betax twiss value at the EXIT of the element.

    ***WARNING!!!: it is always the EXIT, twiss file should be generated accordingly.***
  + How to get betas at any point?

    e.g. `tw.getBeta(3,0.5)` where `3` is the element number and `0.5` is position in the element in meters. `tw.getBeta(3,0)` refers the twiss beta at the ENTRY side.

    To get betas at the output you could use:
```
      tw.elems[3].BETX #(faster option)
      tw.getBeta(3,tw.elems[3].L)
```
  + How to get the phase advance at any point?

    e.g. `tw.getPhase(1,0)`
    exactly the same as in `getBeta()`
### Using fort.18 files or twiss files
  + How do I load a map?

    Create a Map2 from either a fort.18 or a twiss file.
    ``` python
    import mapclass
    mf = mapclass.Map2(4,"fort.18")
    mt = mapclass.Map2(tw,4)
    ```
    where `tw` is a twiss object created with `metaclass2` module and `4` is the map order.

    In the case of a fort.18 file, please, checka that The map order in `Map2` is smaller or equal to the fort.18 map.

  + How to calculate the beamsize

    Using a map `m`, from either a twiss or a fort.18 file (the calculation from a fort.18 is generally more precise)
    ```python
      import math
      var=m.sndmmt('y',[sx,spx,sy,spy,dpp,t])
      mu =m.rstmmt('y',[sx,spx,sy,spy,dpp,t])
      beamsize = math.sqrt(var - mu**2)
    ```
    ***NOTE:***
        `sx`, `spx`, `sy`, `spy` and `t` are one sigma of the gaussian distribution, while
        `dpp`: either one sigma of the gaussian distribution
             or total width of uniform centered energy spread distribution

    ***WARNING:***
        `m.sigma` returns squared value already!
        `m.sigma` and `m.offset` are still valid for backwards compatibility
	
  + How to calculate oide effect?

    `tw.oide(emitn,gamma)`
      where emitn is normalized emittance, gamma is the relativistic factor. 
      Check the function declaration, it has other parameters, this is an example with 
      minimum parameters
  + How to calculate beam size contribution due to radiation in bending
      magnets?

      `tw.sigmaBends2a(E=theenergy)`
      #### NOTE: it returns the squared value.
  + How do I get the horizontal polynomial expression?

    ```python
    print mymap.x   ### prints all the monoms in x
    ```
    It can be used with x, px, y, py, d, s.
  + How do I check an specific component on the map?
    ```python
      print mymap.x[0,0,0,0,1] # 5 dim map, prints the hor. dependence on d
      print mymap.y[0,0,0,1,2] # 5 dim map, prints the  ver. dep. on py * d**2 
    ```
    if the map is 6 dimensional, then, it needs 6 indexes


SIMPLE TROUBLESHOOTING
----------------------
  + 
    ```python
       Traceback (most recent call last):
       File ".../MapClass2/mapclass.py", line 227, in sigma
       sigmaprod = self.__sigma(ind, sig, gaussianDelta)
       File ".../MapClass2/mapclass.py", line 392, in __sigma
       qq = pow(i[5] / dv, ind[5])
    ```
    ***ANSWER:*** Only 5 elements in the passed sigmas list [sx,...], there must be six

+ 
    ```python
    mw=mapclass.Map2("fort.18")
    Traceback (most recent call last):
    File ".../MapClass2/mapclass.py", line 219, in sigma
    for ind1, coeff1 in self[xory].iteritems():
    KeyError: 'x'
    ```
    ***ANSWER:*** Missing argument. A correct call example is
    ` mw=mapclass.Map2(1,"fort.18")`

+ 
    ```python
    print math.sqrt(mw.sigma(xory='x',sig=sigmaFFS))
    Traceback (most recent call last):
    TypeError: can't convert complex to float
    ```
    ***ANSWER:***  When using twiss objects in mapclass, retuned values are complex.

    Use real: `print math.sqrt(mw.sigma(xory='x',sig=sigmaFFS.real)`
+    
    ```python
    No implementation for element:  TWISS MARKER
    ```
    ***ANSWER:*** all is OK, this is an info message


BIBLIOGRAPHY
------------
[1] R. Tomas. Nonlinear optimization of beams lines. PRST - AB 9,081001 (2006).

[2] Martinez, David et al. MAPCLASS2: a code to aid the optimization of 
      lattice design. CERN-ATS-Note-2012-087 TECH. Nov 01, 2012.

[3] Diana Andreea Popescu. Parallel Computing Methods for Particle Accelerator
      Design. Master Thesis. Ecole Polytechnique Federale de Lausanne, CERN.
      July, 2013.

[4] O.R. Blanco et al. CLIC 3TeV Beam Size Optimization with Radiation Effects.
      Proceedings of IPAC2013, Shanghai, China. TUPW0003.

[5] E. FOREST et al. Introduction to Polymorphic Tracking Code. Fibre Bundles,
      Polymorphic Taylor Types and "Exact Tracking". Geneva, Switzerland.
      July 24, 2002.
      CERN-SL-2002-044 (AP).
      KEK-Report 2002-3.

[6] MAD-X. http://mad.web.cern.ch/mad/

[7] R. Tomas. MAPCLASS: a code to optimize high order aberrations.
      CERN-AB-Note-2006-017. Jan 15, 2007.


LICENCE
-------
Code was originally written by R. Tomas (CERN)

CONTACT
--------------------------
email: orblancog@gmail.com, rogelio.tomas@cern.ch
