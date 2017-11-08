MAPCLASS2
=========
### A code to aid the optimisation of lattice design

This Readme.md contains a brief explanation, installation and use instructions.

See, bibliography section :

[1] a conceptual review.

[2] use and detailed python code explanations.

[3] C++ and CUDA libraries.

[4] radiation and Oide effect.

[5] Appendix B, the fort.18 file details.

[6] MAD-X in general

[7] initial MAPCLASS code.

DESCRIPTION
-----------
MAPCLASS2 uses a transfer map either
 * loaded from a fort.18 file generated in MAD-X PTC,
 * or generates a map from a MAD-X twiss table output,

to calculate the first and second order moments of a beam to any given map order. In addition, several map manipulations can be done.

The beam is assumed gaussian in the x, px, y, py and t. The beam energy spread distribution could be chosen to be gaussian or uniform.


INSTALL INSTRUCTIONS
--------------------
Two options are available:
    C++ libraries (recommended) and
    CUDA libraries (some build must be done in the working machine)

1. From a terminal, move to desired installation directory, then do:
```bash
  $ git clone https://github.com/pylhc/MapClass2.git
```
2. Move to the MapClass2 created folder and do:
```bash
  $ git submodule init
  $ git submodule update
```
3. Python variables need to be modified in order to load libraries and run it.
In Linux systems there are three ways (in all cases see below NOTE):
* by adding ~/.bashrc file
```bash
        export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:path_to_libs_boost_1_53_0_libboost_python.so.1.53.0      
        export PYTHONPATH=$PYTHONPATH:path_to_MapClass:path_to_mapbeamline.so
```
* by loading them in the user Python script, for example
```python
      cdll.LoadLibrary("path/to/libs/boost_1_53_0/libboost_python.so.1.53.0")
```
and equivalent for the mapbeamline.so
* by copying the .so files to the working directory
```
    libboost_python.so.1.53.0
    mapbeamline.so
```
#### NOTE: 
*  If you choose C++ LIBRARIES (recommended), use the "C++MapConstruction" dir 
*  If you choose CUDA LIBRARIES (it must be build), use the "CUDAMapConstruction" dir

#### COMMENTS: for the new users of git
To update, move to MapClass2 folder:
```bash
$ git pull
```
To recover in an accidental file remove or modification
```bash
$ git reset --hard
```
SOFTWARE REQUIREMENTS
---------------------
+ Python 2.6 or higher is required (not tested with Python 3.X)
   +  numpy library must be installed or added to the 'libs' folder
+ gcc 4.4.3 or higher.
   +  If working on the cern afs and your gcc is lower than 4.4.3, then do
```bash
      $ source /afs/cern.ch/sw/lcg/external/gcc/4.4/x86_64-slc5/setup.sh
```
+ (OPTIONAL) CUDA libraries

PROVIDED LIBRARIES
------------------
+ boost library for running on lxplus cern server can be found in the libs
    folder or it can be installed when running on local machine


HOW TO USE
----------
From any Python terminal environment or script, it is possible to load the modules mapclass and metaclass2.
```python
>>> import mapclass
>>> import metaclass2
```
now, a twiss or fort.18 file (from MAD-X 5.X.X) is required to perform beamsize
calculations.

### INPUT FILES
+ twiss files from MAD-X
  Twiss file should be generated with at least the following columns:
```
NAME, KEYWORD, S, L, BETX, BETY, ALFX, ALFY, MUX, MUY, DX, DPX, DY, DPY, ANGLE, K1L, K2L, K3L, K4L
```
  ***WARNING: Some functions assume elements referenced to its exit side.***
+ fort.18 files from MAD-X PTC

  Check the PTC module docs in MAD-X as it is very powerful. Here is a common MAD-X coding example to generate fort.18 file.
```
      !!! MAD-X Code
      !PTC  To procude fort.18
        ptc_create_universe;
        ptc_create_layout,model=2,method=6,nst=10;
        ptc_normal,icase=5,no=8,deltap=0.00;
        ptc_end; 
      !!! End MAD-X Code
```

EXAMPLES
--------
Some examples of use can be found in ./doc/FFSexample


QUICK Q/A
---------
### Using twiss files
  + How to load a twiss file?

    In a python script or shell:
```
      import metaclass2
      tw = metaclass2.twiss2("twissfilename")
```
  + How to check what was loaded?

    `len(tw.elems)` tells you how many elements were loaded, while
    `tw.elems[0]` is the first loaded element.

    `tw.elems[0].KEYWORD` is the keyword used in MAD-X for the element zero. In general MAPCLASS2 uses the same 
    twiss file column names except
    for `TNAME` which is `TableNAME` (TNAME does not exist in MAD-X).

    `tw.elems[0].BETX` is the betax twiss value at the EXIT of the element.

    ***WARNING!!!: it is always the EXIT, twiss file should be generated accordingly.***
  + How to get betas at any point?

    e.g. `tw.getBeta(3,0.5)` where `3` is the element number and `0.5` is position in the element. Position equal to 0
    refers to ENTRY side.

    To get betas at the output you could use:
```
      tw.elems[3].BETX #(faster option)
      tw.getBeta(3,tw.elems[3].L)
```
  + How to get the phase advance at any point?

    e.g. `tw.getPhase(#,0)`
    exactly the same as in `getBeta()`
### Using fort.18 files or twiss files
  + How do I load a map?

    Create a Map2 from either a fort.18 or a twiss file.
    ``` python
      import mapclass
      mf = mapclass.Map2(#,"fort.18")
      mt = mapclass.Map2(tw,#)
    ```
    where `tw` is a twiss object created with metaclass2 module
  + How to calculate the beamsize

    Using a map m, either mf or mt (from a fort.18 is more precise)
    ```
    m.sndmmt('y',[sx,spx,sy,spy,dpp,t])
    m.rstmmt('y',[sx,spx,sy,spy,dpp,t])
    beamsize = m.sndmmt - m.rstmmt**2
    ```
    #### NOTE:
        sx,spx,sy,spy and t are one sigma of the gaussian distribution
        dpp: either one sigma of the gaussian distribution
             or total width of uniform centered energy spread distribution
    ####  WARNING:
        `m.sigma` returns squared value already!
        `m.sigma` and `m.offset` are still valid for backwards compatibility
  + How to calculate oide effect?

    `tw.oide(emitn,gamm)`
      where emitn is normalized emittance, gamma is the relativistic factor. 
      Check the function declaration, it has other parameters, this is a 
      minimum implementation
  + How to calculate beam size contribution due to radiation in bending
      magnets?

      `tw.sigmaBends2a(E=theenergy)`
      #### NOTE: it returns the squared value.
  + How do I get the horizontal polynomial expression?

    ```python
    print mymap.x   ### prints all the polynoms in x
    ```
    It can be used with x, px, y, py, d, s.
  + How do I check an specific component on the map?
    ```python
      print mymap.x[0,0,0,0,1] # 5 dim map, prints the hor. dependence on d
      print mymap.y[0,0,0,1,2] # 5 dim map, prints the  ver. dep. on py*d**2 
    ```
    if the map is 6 dimensional, then, it needs 6 indexes


SIMPLE TROUBLESHOOTING
----------------------
### ERRORS
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
    ***ANSWER:*** Missing argument, correct call example is
    ` mw=mapclass.Map2(1,"fort.18")`

+ 
    ```python
    print math.sqrt(mw.sigma(xory='x',sig=sigmaFFS))
    Traceback (most recent call last):
    TypeError: can't convert complex to float
    ```
    ***ANSWER:***  When using twiss objects in mapclass, retuned values are complex.
    Use real:
    ```
    print math.sqrt(mw.sigma(xory='x',sig=sigmaFFS.real)
    ```

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

CONTACT/COMMENTS/QUESTIONS
--------------------------
email: orblancog@gmail.com
