
#include <R.h>
#include <Rmath.h>
#include <Rinternals.h>
#include <Rdefines.h>
#include <R_ext/Applic.h> /* for dgemm */


/**
    number of rows
    *\param x a matrix
*/
int nrow(SEXP x) 
{
    return(INTEGER(getAttrib(x, R_DimSymbol))[0]);
}

/**
    number of columns
    *\param x a matrix
*/
int ncol(SEXP x) 
{
    return(INTEGER(getAttrib(x, R_DimSymbol))[1]);
}

/**
    matrix product x %*% y
    *\param x a matrix
    *\param nrx number of rows of x
    *\param ncx number of cols of x
    *\param y a matrix
    *\param nry number of rows of y
    *\param ncy number of cols of y
    *\param z a matrix of dimension nrx x ncy
*/
void C_matprod(double *x, int nrx, int ncx,
               double *y, int nry, int ncy, double *z)
{
    char *transa = "N", *transb = "N";
    double one = 1.0, zero = 0.0;
    int i;

    if (nrx > 0 && ncx > 0 && nry > 0 && ncy > 0) {
        F77_CALL(dgemm)(transa, transb, &nrx, &ncy, &ncx, &one,
	                x, &nrx, y, &nry, &zero, z, &nrx);
    } else /* zero-extent operations should return zeroes */
	for(i = 0; i < nrx*ncy; i++) z[i] = 0;
}

/**
    trace of the boosting hat operator for each iteration of 
    gradient boosting with componentwise linear base learners 
    (such as smoothing splines) and the final boosting hat operator
    *\param nobs number of observations
    *\param H a list of hat matrices for the base learners
    *\param xselect a vector of selected covariates
*/
SEXP R_trace_gamboost(SEXP nobs, SEXP H, SEXP xselect) 
{
    int n, nn, i, b, B;
    SEXP ans, hatmat, trace;
    double *dhatmat, *dH, *z, *dtrace;
    
    /* allocate memory */
    B = LENGTH(xselect);
    n = INTEGER(nobs)[0];
    nn = n * n;
    z = Calloc(nn, double);
    PROTECT(ans = allocVector(VECSXP, 2));
    SET_VECTOR_ELT(ans, 0, hatmat = allocMatrix(REALSXP, n, n));
    dhatmat = REAL(hatmat);    
    for (i = 0; i < nn; i++) dhatmat[i] = 0.0;
    SET_VECTOR_ELT(ans, 1, trace = allocVector(REALSXP, B));
    dtrace = REAL(trace);

    /* for each boosting iteration */    
    for (b = 0; b < B; b++) {

        /* hat matrix of base learner */
        dH = REAL(VECTOR_ELT(H, INTEGER(xselect)[b] - 1));
        
        /* update boosting hat matrix */
        C_matprod(dH, n, n, dhatmat, n, n, z);
        for (i = 0; i < nn; i++) 
            dhatmat[i] += dH[i] - z[i];

        /* its trace */
        dtrace[b] = 0.0;
        for (i = 0; i < n; i++) 
            dtrace[b] += dhatmat[i + n * i];
            
    }
    Free(z);
    UNPROTECT(1);
    return(ans);
}

/* z := t(x) %*% (I - B) */
void C_dvecImat(double *x, int n, int xselect, double *B, double *z)
{

    int i, j, jn, xn = xselect * n;
    
    for (j = 0; j < n; j++) {
        z[j] = 0.0;
        jn = j * n;

        for (i = 0; i < j; i++)
            z[j] -= x[xn + i] * B[jn + i];

        /* diagonal elements */
        i = j;
        z[j] += x[xn + i] * (1 - B[jn + i]);

        for (i = (j + 1); i < n; i++)
            z[j] -= x[xn + i] * B[jn + i];
    }
}
            
/* B := B + xf %*% z */
void C_updateB(double *B, int n, double *xf, int xselect, double *z)
{

    int i, j, jn, xn = xselect * n;
    
    for (j = 0; j < n; j++) {
        jn = j * n;
        for (i = 0; i < n; i++)
            B[jn + i] += xf[xn + i] * z[j];
    }
}

/**
    trace of the boosting hat operator for each iteration of
    gradient boosting with componentwise linear models and
    the final boosting hat operator
    *\param x design matrix
    *\param MPinv t( nu * (t(x) %*% x)^{-1} x
    *\param xselect B-vector of selected covariates
*/
SEXP R_trace_glmboost(SEXP x, SEXP MPinv, SEXP xselect) {

    SEXP ans, hatmatrix, trace;
    int i, B, b, n, *ixselect;
    double *dhatmatrix, *dtrace, *z;
    
    B = LENGTH(xselect);
    n = nrow(x);
    
    /* allocate memory */
    ixselect = INTEGER(xselect);
    PROTECT(ans = allocVector(VECSXP, 2));
    SET_VECTOR_ELT(ans, 0, hatmatrix = allocMatrix(REALSXP, n, n));
    dhatmatrix = REAL(hatmatrix);
    SET_VECTOR_ELT(ans, 1, trace = allocVector(REALSXP, B));
    dtrace = REAL(trace);
    for (i = 0; i < n * n; i++) dhatmatrix[i] = 0.0;
    z = Calloc(n, double);
    
    /* for each boosting iteration */
    for (b = 0; b < B; b++) {

        /* update boosting hat operator */
        C_dvecImat(REAL(MPinv), n, ixselect[b] - 1, dhatmatrix, z);
        C_updateB(dhatmatrix, n, REAL(x), ixselect[b] - 1, z);
        
        /* its trace for the bth iteration */
        dtrace[b] = 0.0;
        for (i = 0; i < n; i++)
            dtrace[b] += dhatmatrix[i + i * n];
        
    }
    UNPROTECT(1);
    Free(z);
    return(ans);
}

/**
    trace of a quadratic matrix
    *\param B quadratic matrix
*/

SEXP sumdiag(SEXP B) {

    SEXP ans;
    double *dans, *dB;
    int i, n;
    
    PROTECT(ans = allocVector(REALSXP, 1));
    dans = REAL(ans);
    dans[0] = 0;
    dB = REAL(B);
    n = nrow(B);
    
    for (i = 0; i < n; i++) 
        dans[0] += dB[i * n + i];
    
    UNPROTECT(1);
    return(ans);
}

/**
    negative gradient of the partial likelihood of a Cox model
    see formula (4.65) in Chapter 4 of Greg Ridgeway's thesis    
    http://www.i-pensieri.com/gregr/papers/thesis.pdf
    and Section 4.5 in vignette("gbm")
    *\param time survival times
    *\param event censoring indicate (event == 1 means dead)
    *\param f boosting fit
    *\param w weights
*/
SEXP ngradientCoxPLik(SEXP time, SEXP event, SEXP f, SEXP w) {

    SEXP ans;
    double *dtime, *df, *dans, *dummy, *dw;
    int *ievent, i, j, k, n;
    
    /* we don't assume the variables to be ordered w.r.t. time */

    /* allocate memory */
    n = LENGTH(time);
    PROTECT(ans = allocVector(REALSXP, n));
    dans = REAL(ans);
    dtime = REAL(time);
    ievent = INTEGER(event);
    df = REAL(f);
    dw = REAL(w);
    dummy = Calloc(n, double);
    
    for (i = 0; i < n; i++) {
        df[i] = exp(df[i]);
        dans[i] = 0.0;
    }
        
    for (j = 0; j < n; j++) {
        for (k = 0; k < n; k++) {
            if (!(dtime[j] > dtime[k]) || j == k)
                dummy[j] += dw[k] * df[k];
        }
    }
    
    for (i = 0; i < n; i++) {
        for (j = 0; j < n; j++) {
            if (ievent[j] & (!(dtime[j] > dtime[i]))) 
                dans[i] += dw[j] * df[i] / ((dummy[j] == 0) ? 1 : dummy[j]);
        }
        dans[i] = ievent[i] - dans[i];
    }

    Free(dummy);
    UNPROTECT(1);
    return(ans);
}

/**
    copy an R object
    *\param x an R object
*/
SEXP copymem (SEXP x) {
    
    return(duplicate(x));
}


int C_intmax(const int *x, const int n) {
   int tmp = 0;
   int i;
   
   for (i = 0; i < n; i++) {
       if (x[i] == NA_INTEGER) continue;
       if (x[i] > tmp) tmp = x[i];
   }
   return(tmp);
}


/**
    tapply(y, index, sum)
*/
SEXP R_ysum (SEXP y, SEXP index) {

    int n, m, i;
    SEXP ans;
    double *dans, *dy;
    int *iindex;
    
    n = LENGTH(y);
    dy = REAL(y);
    if (n != LENGTH(index)) error("dimensions don't match");
    iindex = INTEGER(index);
    m = C_intmax(iindex, n);

    PROTECT(ans = allocVector(REALSXP, m));
    dans = REAL(ans);
    for (i = 0; i < m; i++) {
        dans[i] = 0.0;
    }

    for (i = 0; i < n; i++) {
        if (iindex[i] == NA_INTEGER) continue;
        dans[iindex[i] - 1] += dy[i];
    }
    UNPROTECT(1);
    return(ans);
}


/**
    partial likelihood of a Cox model
    see formula (4.62) in Chapter 4 of Greg Ridgeway's thesis
    http://www.i-pensieri.com/gregr/papers/thesis.pdf
    and Section 4.5 in vignette("gbm")
    *\param time survival times (ordered!)
    *\param expf exp(boosting fit)
*/

SEXP R_risk (SEXP time, SEXP expf) {

    SEXP ans;
    double *dtime, *dexpf, *dans;
    int i, j, n;
    
    n = LENGTH(time);
    PROTECT(ans = allocVector(REALSXP, n));
    dans = REAL(ans);
    dtime = REAL(time);
    dexpf = REAL(expf);
    
    for (i = 0; i < n; i++) {
        dans[i] = 0;
        for (j = 0; j < n; j++)
            if (!(dtime[j] < dtime[i]) || j == i) 
                dans[i] += dexpf[j];
    }
    
    UNPROTECT(1);
    return(ans);
}

SEXP R_mcumsum (SEXP x) {

    SEXP ans;
    double *dans, *dx;
    int nr, nc, i, j, jj;
    
    nr = nrow(x);
    nc = ncol(x);

    PROTECT(ans = allocMatrix(REALSXP, nr, nc));
    dans = REAL(ans);
    dx = REAL(x);

    for (j = 0; j < nc; j++) {
        for (i = 0; i < nr; i++) {
            dans[j * nr + i] = 0.0;
        }
    }

    for (jj = 0; jj < nc; jj++) {
        for (j = jj; j < nc; j++) {
            for (i = 0; i < nr; i++) {
                dans[j * nr + i] += dx[jj * nr + i];
            }
        }
    }
    UNPROTECT(1);
    return(ans);
}
