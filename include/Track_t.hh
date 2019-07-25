#ifndef TRACK_T_HH_
#define TRACK_T_HH_ 1


// #include
#include <stdlib.h>

// #include <vector>
// #include "TROOT.h"
// gROOT->ProcessLine("#include <vector>");

// Track structure
typedef struct {
  double pX;
  double pY;
  double pZ;
  int    ID;
  int   PDG;
}Track_t;

#endif /* TRACK_T_HH_ */
