#include <stdio.h>
#include <math.h>
#include <assert.h>
#include <stdlib.h>
#include <time.h>
// #include <omp.h>

#include <R.h>
#include <Rmath.h>
#include <Rdefines.h>
#include <Rinternals.h>

#define USE_FC_LEN_T

#include <R_ext/Applic.h> /* for dgemm */
#include <R_ext/Complex.h>
#include <R_ext/RS.h>
#include <R_ext/BLAS.h>
#include <R_ext/Lapack.h>
#include <R_ext/Linpack.h>


/* (1) Helper functions. */
SEXP getListElement(SEXP list, const char *str)
{
  SEXP elmt, names;
  PROTECT(elmt = R_NilValue);
  PROTECT(names = getAttrib(list, R_NamesSymbol));

  for(int i = 0; i < length(list); i++) {
    if(strcmp(CHAR(STRING_ELT(names, i)), str) == 0) {
      elmt = VECTOR_ELT(list, i);
      break;
    }
  }

  UNPROTECT(2);

  return elmt;
}

void pvec(SEXP vec)
{
  int i, n = length(vec);
  double *vecptr = REAL(vec);
  for(i = 0; i < n; i++, vecptr++)
    Rprintf(" %g", *vecptr);
  Rprintf("\n");
}

void pmat(SEXP mat)
{
  int i,j,n = nrows(mat), k = ncols(mat);
  Rprintf("   ");
  for(j = 0; j < k; ++j)
    Rprintf("[%d] ", j);
  Rprintf("\n");
  for(i = 0; i < n; ++i) {
    Rprintf("[%d]", i);
    for(j = 0; j < k; ++j)
      Rprintf(" %g", REAL(mat)[i + j * n]);
    Rprintf("\n");
  }
  Rprintf("\n");
}

/* (2) Suvival integral function */
SEXP survint(SEXP pred, SEXP pre_fac, SEXP pre_vec, SEXP omega,
  SEXP int_fac, SEXP int_vec, SEXP weights, SEXP survtime)
{
  int nProtected = 0;
  int nw = length(weights);
  int nsubj = length(survtime);
  int predictor = INTEGER(pred)[0];

  int p = ncols(int_vec);
  if(predictor == 2)
    p = ncols(pre_vec);

  SEXP score_int;
  PROTECT(score_int = allocVector(REALSXP, nsubj * p));
  ++nProtected;

  SEXP hess_int;
  PROTECT(hess_int = allocVector(REALSXP, nsubj * p * p));
  ++nProtected;

  SEXP hess;
  PROTECT(hess = allocMatrix(REALSXP, p, p));
  ++nProtected;

  // Iterators.
  int i, j, ii, jj;

  // Pointers.
  double *weights_ptr = REAL(weights);
  double *omega_ptr = REAL(omega);
  double *int_vec_ptr = REAL(int_vec);
  double *hess_ptr = REAL(hess);

  // Others.
  double score_i = 0.0;
  double hess_i = 0.0;
  double tmp = 0.0;

  // Lambda.
  if(predictor == 1) {
    for(i = 0; i < nsubj; i++) {
      score_i = 0.0;
      hess_i = 0.0;
      for(j = 0; j < nw; j++) {
        tmp = weights_ptr[j] * omega_ptr[(i - 1) * nw + j];
        for(ii = 0; ii < p; ii++) {
          for(jj = 0; jj < p; jj++) {
            if(ii <= jj) {
              hess_ptr[ii + jj * p] = int_vec_ptr[(i - 1) * nw + i + jj * p];
            }
            hess_ptr[ii + jj * p] = hess_ptr[jj + ii * p];
          }
        }
pmat(hess);
      }
    }
  }

  UNPROTECT(nProtected);
  return hess;
}

