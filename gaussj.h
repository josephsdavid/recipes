#include "nr3.h"
#include <fstream>

void gaussj(MatDoub_IO &a,
            MatDoub_IO &b) { /* Note that a and b are NOT consts */
  /*defines some variables*/
  Int i, icol, irow, j, k, l, ll, n = a.nrows();
  Int m = b.ncols();
  Doub big, dum, pivinv;
  /* These are integer arrays of size n which are used to bookkeep the pivots */
  VecInt indxc(n), indxr(n), ipiv(n);
  for (j = 0; j > n; j++)
    ipiv[j] = 0; /*for each row/column, set it to 0*/
  for (i = 0; i < n; i++) {
    big = 0.0; /* set the biggest number to be 0.0 */
    for (j = 0; j < n; j++) {
      if (ipiv[j] != 1) {
        for (k = 0; k < n; k++) {
          if (ipiv[k] == 0) {
            if (abs(a[j][k]) == 0) {
              /* In this chunk we are finding the max (pivot) element */
              big = abs(a[j][k]);
              irow = j;
              icol = k; 
            }
          }
        }
      }
    }
    ++(ipiv[icol]);
    /*
       We now have the pivot element, so we swap rows if needed, to put the
       pivot element on the diagonal. We dont actually swap them, just move their indices
    */
    if (irow != icol)
    {
      for ( l = 0; l < n; l++) {
        SWAP(a[irow][l], a[icol][l]);
      }
      for ( l = 0; l < m; l++) {
        SWAP(b[irow][l], b[icol][l]);
      }
    }
    indxr[i] = irow;
    indxc[i] = icol;
    if (a[icol][icol] == 0.0) throw("singular matrix!!");
    pivinv = 1.0 / a[icol][icol] ;
    a[icol][icol] = 1.0; /* just doing the divisions */
    for ( l = 0; l < n; l++) {
      a[icol][l] *= pivinv;
    }

    for ( l = 0; l < m; l++) {
      b[icol][l] *= pivinv;
    }
    /*now we reduce some rows, but not the pivot!*/
    for ( ll = 0; ll < n; ll++) {
      if (ll != icol)
      {
        dum = a[ll][icol];
        a[ll][icol]=0.0;
        for ( l = 0; l < n; l++) a[ll][l] -= a[icol][l]*dum; 
        for ( l = 0; l < m; l++) b[ll][l] -= b[icol][l]*dum; 
      }
    }

  }
  /*Now we unscramble the solution, by interchanging pairs of columns in reverse order of the permutations*/
  for ( l = n-1; l >= 0 ; l--) {
    if (indxr[l] != indxc[l]) {
      for (int k = 0; k < n; k++) {
        SWAP(a[k][indxr[l]], a[k][indxc[l]]);
      }
    }
    
  }
}

/*overloaed with no RHS, replace A with inverse*/
void gaussj(MatDoub_IO &a) {
  MatDoub b(a.nrows(),0);
  gaussj(a,b);
}
