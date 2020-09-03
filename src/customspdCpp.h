#ifndef CSDPDECLARATIONS
#define CSDPDECLARATIONS

#include "declarations.h"

extern "C" int custom_sdpCpp(
     int n,
     int k,
     struct blockmatrix C,
     double *a,
     struct constraintmatrix *constraints,
     double constant_offset,
     struct blockmatrix *pX,
     double **py,
     struct blockmatrix *pZ,
     double *ppobj,
     double *pdobj);

#endif // CSDPDECLARATIONS

