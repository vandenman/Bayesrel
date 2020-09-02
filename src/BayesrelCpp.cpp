
//#include <R.h>
//#include <Rinternals.h>
//#include <stdlib.h> // for NULL
//#include <R_ext/Rdynload.h>
#include <declarations.h>
#include <stdio.h>
#include <RcppArmadillo.h>
#include <Rcpp.h>
#include "customsdp.h"


//[[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;
using namespace arma;

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

//[[Rcpp::export]]
Rcpp::List csdpArma(
              int n_p,
              int nconstraints_p,
              int nblocks_p,
              arma::ivec blocktypes_p,
              arma::ivec blocksizes_p,
              Rcpp::List C_p,
              Rcpp::List A_p,
              arma::dvec b_p)
{


    struct blockmatrix C;
    struct constraintmatrix *constraints;
    struct blockmatrix X, Z;
    double *y, *b;
    double pobj, dobj;
    int status;


    /*
     * setup C
     */
    C = blkmatrix_R2csdpArma(C_p);

    /*
     * setup constraints
     */
    constraints = constraints_R2csdpArma(A_p);

    /*
     * Allocate storage for RHS
     */
    b = double_vector_R2csdpArma(nconstraints_p,b_p);

    /*
     * Create an initial solution. This allocates space for X, y, and Z,
     * and sets initial values
     */
    initsoln(n_p,nconstraints_p,C,b,constraints,&X,&y,&Z);

    /*
     * Solve the problem
     */
    status = custom_sdpCpp(n_p,nconstraints_p,C,b,constraints,0.0,&X,&y,&Z,&pobj,&dobj);

    /*
     * Grab the results
     */

    /*
     * Grab X
     */
    Rcpp::List X_p = blkmatrix_csdp2RArma(X);

    /*
     * Grab Z
     */
    Rcpp::List Z_p = blkmatrix_csdp2RArma(Z);

    /* Copy y */
    arma::dvec y_p = double_vector_csdp2RArma(nconstraints_p, y);


    free_prob(n_p,nconstraints_p,C,b,constraints,X,y,Z);

    Rcpp::List ret;
    ret.push_back(X_p);
    ret.push_back(Z_p);
    ret.push_back(y_p);
    ret.push_back(pobj);
    ret.push_back(dobj);
    ret.push_back(status);

  return ret;
}


//static const R_CallMethodDef CallEntries[] = {
//  {"csdpArma",          (DL_FUNC) &csdpArma,          8},
//  {NULL, NULL, 0}
//};
//
//void R_init_Bayesrel(DllInfo *dll)
//{
//  R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
//  R_useDynamicSymbols(dll, FALSE);
//}
