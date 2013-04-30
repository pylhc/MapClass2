setenv gcc_config_version 4.4.3
setenv LCGPLAT x86_64-slc5-gcc44-opt
setenv LCG_lib_name lib64

setenv LCG_contdir /afs/cern.ch/sw/lcg/contrib
setenv LCG_gcc_home ${LCG_contdir}/gcc/${gcc_config_version}/${LCGPLAT}

setenv PATH ${LCG_gcc_home}/bin:${PATH}


setenv LD_LIBRARY_PATH ${LCG_gcc_home}/${LCG_lib_name}:${LD_LIBRARY_PATH}


