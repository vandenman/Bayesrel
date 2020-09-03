//
//  customsdp.h
//
//
//  Created by Julius Pfadt on 02.09.20.
//

#ifndef customsdp_h
#define customsdp_h
#include <RcppArmadillo.h>


arma::ivec int_vector_csdp2RArma(int n, int *y);

arma::dvec double_vector_csdp2RArma(int n, double *y);

int * int_vector_R2csdpArma(int n, arma::ivec y);

double * double_vector_R2csdpArma(int n, arma::dvec y);

struct blockmatrix blkmatrix_R2csdpArma(Rcpp::List X);

Rcpp::List blkmatrix_csdp2RArma(struct blockmatrix X);

struct constraintmatrix *constraints_R2csdpArma(Rcpp::List A);

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


#endif /* customsdp_h */
