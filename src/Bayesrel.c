
#include <R.h>
#include <Rinternals.h>
#include <declarations.h>
#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>

SEXP double_vector_csdp2R(int, double*);
double *double_vector_R2csdp(int, SEXP);
struct blockmatrix blkmatrix_R2csdp(SEXP);
SEXP blkmatrix_csdp2R(struct blockmatrix);
struct constraintmatrix *constraints_R2csdp(SEXP);
SEXP int_vector_csdp2R(int, int*);
SEXP double_vector_csdp2R(int, double*);
int *int_vector_R2csdp(int, SEXP);

int custom_sdp(
    int,
    int,
    struct blockmatrix,
    double *,
    struct constraintmatrix *,
    double,
    struct blockmatrix *,
    double **,
    struct blockmatrix *,
    double *,
    double *);
/*
 int custom_sdpCpp(
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
 */

SEXP csdp(SEXP n_p,
	  SEXP nconstraints_p,
	  SEXP nblocks_p,
	  SEXP blocktypes_p,
	  SEXP blocksizes_p,
	  SEXP C_p,
	  SEXP A_p,
	  SEXP b_p)
{
  SEXP X_p, Z_p, y_p, ret, pobj_p, dobj_p, status_p;
  enum AIJ_SLOTS {AIJ_NNZ, AIJ_IIND, AIJ_JIND, AIJ_X};

  int n, nblocks, nconstraints, *blocktypes, *blocksizes;

  struct blockmatrix C;
  struct constraintmatrix *constraints;
  struct blockmatrix X, Z;
  double *y, *b;
  double pobj, dobj;
  int status;

  n = INTEGER(n_p)[0];
  nblocks = INTEGER(nblocks_p)[0];
  nconstraints = INTEGER(nconstraints_p)[0];
  blocktypes = INTEGER(blocktypes_p);
  blocksizes = INTEGER(blocksizes_p);

  /*
   * setup C
   */
  C = blkmatrix_R2csdp(C_p);

  /*
   * setup constraints
   */
  constraints = constraints_R2csdp(A_p);

  /*
   * Allocate storage for RHS
   */
  b = double_vector_R2csdp(nconstraints,b_p);
  if (b==NULL) error("Failed to allocate storage for RHS vector b!\n");

  /*
   * Create an initial solution. This allocates space for X, y, and Z,
   * and sets initial values
   */
  initsoln(n,nconstraints,C,b,constraints,&X,&y,&Z);

  /*
   * Solve the problem
   */
  status = custom_sdp(n,nconstraints,C,b,constraints,0.0,&X,&y,&Z,&pobj,&dobj);

  /*
   * Grab the results
   */

  /*
   * Grab X
   */
  X_p = PROTECT(blkmatrix_csdp2R(X));

  /*
   * Grab Z
   */
  Z_p = PROTECT(blkmatrix_csdp2R(Z));

  /* Copy y */
  y_p = PROTECT(double_vector_csdp2R(nconstraints, y));

  pobj_p = PROTECT(allocVector(REALSXP,1)); REAL(pobj_p)[0] = pobj;
  dobj_p = PROTECT(allocVector(REALSXP,1)); REAL(dobj_p)[0] = dobj;
  status_p = PROTECT(allocVector(INTSXP,1)); INTEGER(status_p)[0] = status;

  free_prob(n,nconstraints,C,b,constraints,X,y,Z);

  ret = PROTECT(allocVector(VECSXP,6));
  SET_VECTOR_ELT(ret,0,X_p);
  SET_VECTOR_ELT(ret,1,Z_p);
  SET_VECTOR_ELT(ret,2,y_p);
  SET_VECTOR_ELT(ret,3,pobj_p);
  SET_VECTOR_ELT(ret,4,dobj_p);
  SET_VECTOR_ELT(ret,5,status_p);

  UNPROTECT(7);
  return ret;
}


static const R_CallMethodDef CallEntries[] = {
  {"csdp",          (DL_FUNC) &csdp,          8},
  {NULL, NULL, 0}
};

void R_init_Bayesrel(DllInfo *dll)
{
  R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
  R_useDynamicSymbols(dll, FALSE);
}
