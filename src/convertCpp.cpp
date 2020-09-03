//
//
//
//  Created by Julius Pfadt on 26.08.20.
//

#include <stdio.h>
#include <RcppArmadillo.h>
#include <Rcpp.h>
extern "C" {
#include "blockmat.h"
}


//[[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;
using namespace arma;

ivec int_vector_csdp2RArma(int n, int *y)
{
    return ivec(y, n+1);
}

dvec double_vector_csdp2RArma(int n, double *y)
{
    return dvec(y, n+1);
}

int * int_vector_R2csdpArma(int n, ivec y)
{
  // return y.memptr(); // <- this should also work
  int *ret;
  int i;
  ret = (int *) malloc((n+1) * sizeof(int));
  if (ret == NULL)
    return NULL;
  for (i=1; i<=n; i++)
    ret[i] = y[i];
  return ret;
}



double * double_vector_R2csdpArma(int n, dvec y)
{
  // return y.memptr(); // <- this should also work
  double *ret;
  int i;
  ret = (double *) malloc((n+1) * sizeof(double));
  if (ret == NULL)
    return NULL;
  for (i=1; i<=n; i++)
    ret[i] = y[i];
  return ret;
}


/*
 the input of this function is a kind of list but of class csdpBlkMat
 dont know if a Rcpp type list can also work
*/
struct blockmatrix blkmatrix_R2csdpArma(List X)
{
    struct blockmatrix ret;
    int nblocks = X["nblocks"];
    List blocks = X["blocks"];
    ret.nblocks = nblocks;
    ret.blocks = (struct blockrec *) malloc((nblocks + 1) * sizeof(struct blockrec));
    for (int j=1; j<=nblocks; j++) {
      List cur_block = blocks[j-1];
      int blksize = cur_block["blocksize"];
      ret.blocks[j].blocksize = blksize;
      int blktype = cur_block["blockcategory"];
      ret.blocks[j].blockcategory = (blktype == 1) ? MATRIX : DIAG;
      if (blktype == 1) {
        int allocsize = blksize*blksize;
        ret.blocks[j].data.mat = (double *) malloc(allocsize * sizeof(double));
        dvec dblvec = cur_block["data"];
        for (int k=0; k<allocsize; k++)
            ret.blocks[j].data.mat[k] = dblvec[k];
      }
      else {
        ret.blocks[j].data.vec = double_vector_R2csdpArma(blksize, cur_block["data"]);
      }
    }
    return ret;
}


List blkmatrix_csdp2RArma(struct blockmatrix X)
{
    List ret;
    dvec data;
    List blocks;

    int j,k, allocsize;

    int nblocks = X.nblocks;
    ret.push_back(nblocks, "n.blocks");

    for (j=1; j<=X.nblocks; j++) {
        int blocksize = X.blocks[j].blocksize;
        int blockcategory = (X.blocks[j].blockcategory == MATRIX) ? 1 : 2;

        if (X.blocks[j].blockcategory == MATRIX) {
            allocsize = X.blocks[j].blocksize * X.blocks[j].blocksize;
            for (k=0; k<allocsize; k++)
                data[k] = X.blocks[j].data.mat[k];
        } else {
            data = double_vector_csdp2RArma(X.blocks[j].blocksize, X.blocks[j].data.vec);
        }
        blocks.insert(j-1, List::create(Named("blocksize")=blocksize, _["blockcategory"]=blockcategory,
                                        _["data"] = data));

    }
    ret.push_back(blocks, "blocks");
    return ret;
}


struct constraintmatrix *constraints_R2csdpArma(List A)
{
    struct constraintmatrix *constraints;
    struct sparseblock *blockptr;

    int nconstraints, nblocks;
    List Ai, Aij;

    int i,j;

    nconstraints = A.length();
    constraints = (struct constraintmatrix *) malloc((nconstraints + 1) * sizeof(struct constraintmatrix));

    for (i=1; i<=nconstraints; i++) {
        Ai = A[i-1];
        /*
         * Terminate block linked list with NULL
         */
        constraints[i].blocks = NULL;
        nblocks = Ai.length();

        for (j=nblocks; j>=1; j--) {
            Aij = Ai[j-1];
            /*
             * Allocate block data structure
             */
            blockptr = (struct sparseblock *) malloc(sizeof(struct sparseblock));
            /*
             * Initialize block data structure
             */
            blockptr->blocknum=Aij["blocknum"];
            blockptr->blocksize=Aij["blocksize"];
            blockptr->constraintnum=Aij["constraintnum"];
            blockptr->next=NULL;
            blockptr->nextbyblock=NULL;
            blockptr->numentries=Aij["numentries"];

            /*
             * Enter data
             */
            blockptr->iindices = int_vector_R2csdpArma(blockptr->numentries, Aij["iindices"]);

            blockptr->jindices = int_vector_R2csdpArma(blockptr->numentries, Aij["jindices"]);

            blockptr->entries = double_vector_R2csdpArma(blockptr->numentries, Aij["entries"]);

            /*
             * Insert block into linked list
             */
            blockptr->next=constraints[i].blocks;
            constraints[i].blocks=blockptr;
        }
    }
  return constraints;
}

