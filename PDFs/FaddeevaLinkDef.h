//#ifndef EXTLIBDMY_H
//#define EXTLIBDMY_H

#ifdef __CINT__

#pragma link off all globals;
#pragma link off all classes;
#pragma link off all functions;
//#ifndef __MAKECINT__
// Include extlib's original header file here.
// This part is only read by C++ compiler.
//#include "Faddeeva.hh"
//#else  // __MAKECINT__
// Declare interface to the external library.
// This part is read by cint/makecint only and
// must be described within cint limitation.
// #pragma link statements can be added as option.
//class ExtLib { ... };
//#pragma link C++ class ExtLib;
//#pragma link C++ global Faddeeva;
#pragma link C++ function erfc;
//#endif  // __MAKECINT__
#endif


//#endif // EXTLIBDMY_H
