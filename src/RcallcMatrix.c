#include <R.h>
#include <Rinternals.h>
#include <Rdefines.h>
#include <R_ext/Utils.h>
#include <R_ext/BLAS.h>


#ifdef USING_R
#define ALLOC(a,b)  R_alloc(a,b)
#else
#define ALLOC(a,b)  S_alloc(a,b)
#endif

void doublematrix(double *array, int ncol, int nrow, double ** pointer);
int cholesky3(double ** matrix, int n, double *diag, double toler,double * y);
int cholesky3R(SEXP matrix,double *diag, double toler,SEXP y);


SEXP quicksolve1(SEXP mat0,SEXP Yresp) {/*m(row)*n(col)*/
  SEXP beta,mat1;
  SEXP out_list,out_list_names;
  SEXP S_flag;
  double ** Mat1;
  double * beta_i;
  int n,m;
  int i=0,j=0;
  int flag=-1;
  
  __PRE:
  PROTECT(mat1=duplicate(mat0));
  m=INTEGER(getAttrib(mat1,R_DimSymbol))[0];
  n=INTEGER(getAttrib(mat1,R_DimSymbol))[1];
  /*PROTECT( out_mat =allocMatrix(REALSXP,m,n) );*/
  Mat1 = (double **) ALLOC(m, sizeof(double *));
  doublematrix(REAL(mat1),n,m,Mat1);
  PROTECT( beta=duplicate(Yresp) );
  beta_i=REAL(beta);
  PROTECT(S_flag   = allocVector(INTSXP, 1));
  
  /*flag=cholesky3R(mat1, NULL , 1e-10,beta);*/
  flag=cholesky3((double **)Mat1,(int)n,(double *) NULL , 1e-9,beta_i);
  INTEGER(S_flag)[0]=(int) flag;	  
  __DONE:
  /*
  construct return list
  */
  PROTECT(out_list=allocVector(VECSXP,2));
  SET_VECTOR_ELT(out_list, 0, beta);
  SET_VECTOR_ELT(out_list, 1, S_flag);
   
  PROTECT(out_list_names = allocVector(STRSXP, 2));
    SET_STRING_ELT(out_list_names, 0, mkChar("beta"));
    SET_STRING_ELT(out_list_names, 1, mkChar("Convergences"));
    
#ifdef USING_R
    setAttrib(out_list, R_NamesSymbol, out_list_names);
#else
    SET_NAMES(out_list, out_list_names);
#endif      
  UNPROTECT(5);/*matrix+beta+ 2 list*/
  return(out_list);
}

void doublematrix(double *array, int ncol, int nrow, double ** pointer){
    int i;
    for (i=0; i<nrow; i++) {
	pointer[i] = array;
	array += ncol;
	}
/*	printf("%4f\t%4f\n",pointer[nrow-1][nrow-2],pointer[nrow-1][nrow-1]);*/
    }


int cholesky3(double ** matrix, int n, double *diag, double toler,double * y){

    double temp;
    int  i,j,k;
    double eps, pivot;
    int rank;
    int n2;
    int nonneg;
	int m=0;
	
	
    n2 = n-m;    /* number of full covariates */
       
    nonneg=1;
    eps =0;
    for (i=0; i<m; i++) if (diag[i] <eps) eps = diag[i];
    for (i=0; i<n2; i++) if (matrix[i][i+m] > eps)  eps = matrix[i][(i+m)]; 
    eps *= toler;


    rank =0;
    /* pivot out the diagonal elements */
    for (i=0; i<m; i++) {
	pivot = diag[i];
        if (pivot < eps) {
            for (j=0; j<n2; j++) matrix[j,i] =0;
            if (pivot < -8*eps) nonneg= -1;
            }
	else {
	    rank++;
	    for (j=0; j<n2; j++) {
		temp = matrix[j][i] / pivot;
		matrix[j][i] = temp;
		matrix[j][(j+m)] -= temp*temp*pivot;
		for (k=(j+1); k<n2; k++) matrix[k][(j+m)] -= temp*matrix[k][i];
		}
	    }
        }

    /* Now the rest of the matrix */
    for (i=0; i<n2; i++) {
	pivot = matrix[i][(i+m)];
	if (pivot < eps) {
	    for (j=i; j<n2; j++) matrix[j][(i+m)] =0;  /* zero the column */
            if (pivot < -8*eps) nonneg= -1;
	    }
	else  {
	    rank++;
	    for (j=(i+1); j<n2; j++) {
		temp = matrix[j][(i+m)]/pivot;
		matrix[j][(i+m)] = temp;
		matrix[j][(j+m)] -= temp*temp*pivot;
		for (k=(j+1); k<n2; k++) matrix[k][(j+m)] -= temp*matrix[k][(i+m)];
		}
	    }
	}
	
	CHOLSKY2:
     /*
     ** solve Fb =y
     */
     for (i=0; i<n; i++) {
	  temp = y[i] ;
	  for (j=0; j<i; j++)
	       temp -= y[j] * matrix[i][j] ;
	  y[i] = temp ;
	  }
     /*
     ** solve DF'z =b
     */
     for (i=(n-1); i>=0; i--) {
	  if (matrix[i][i]==0)  y[i] =0;
	  else {
	      temp = y[i]/matrix[i][i];
	      for (j= i+1; j<n; j++)
		   temp -= y[j]*matrix[j][i];
	      y[i] = temp;
	      }
	  }
	
    return(rank * nonneg);
    }

int cholesky3R(SEXP matrix, double *diag, double toler,SEXP y){

    double temp;
    int  i,j,k;
    double eps, pivot;
    int rank;
    int n2;
    int nonneg;
	int n; 
	int m=0;
	
	n=INTEGER(getAttrib(matrix,R_DimSymbol))[0];
	
    n2 = n-m;    /* number of full covariates */
   
    nonneg=1;
    eps =0;
    for (i=0; i<m; i++) if (diag[i] <eps) eps = diag[i];
    for (i=0; i<n2; i++) if (REAL(matrix)[i+(i+m)*n] > eps)  eps = REAL(matrix)[i+(i+m)*n];
    eps *= toler;

    rank =0;
    /* pivot out the diagonal elements */
    for (i=0; i<m; i++) {
	pivot = diag[i];
        if (pivot < eps) {
            for (j=0; j<n2; j++) REAL(matrix)[j+i*n] =0;
            if (pivot < -8*eps) nonneg= -1;
            }
	else {
	    rank++;
	    for (j=0; j<n2; j++) {
		temp = REAL(matrix)[j+i*n] / pivot;
		REAL(matrix)[j+i*n] = temp;
		REAL(matrix)[j+(j+m)*n] -= temp*temp*pivot;
		for (k=(j+1); k<n2; k++) REAL(matrix)[k+(j+m)*n] -= temp*REAL(matrix)[k+i*n];
		}
	    }
        }

    /* Now the rest of the matrix */
    for (i=0; i<n2; i++) {
	pivot = REAL(matrix)[i+(i+m)*n];
	if (pivot < eps) {
	    for (j=i; j<n2; j++) REAL(matrix)[j+(i+m)*n] =0;  /* zero the column */
            if (pivot < -8*eps) nonneg= -1;
	    }
	else  {
	    rank++;
	    for (j=(i+1); j<n2; j++) {
		temp = REAL(matrix)[j+(i+m)*n]/pivot;
		REAL(matrix)[j+(i+m)*n] = temp;
		REAL(matrix)[j+(j+m)*n] -= temp*temp*pivot;
		for (k=(j+1); k<n2; k++) REAL(matrix)[k+(j+m)*n] -= temp*REAL(matrix)[k+(i+m)*n];
		}
	    }
	}
	
	CHOLSKY2:
     /*
     ** solve Fb =y
     */
     for (i=0; i<n; i++) {
	  temp = REAL(y)[i] ;
	  for (j=0; j<i; j++)
	       temp -= REAL(y)[j] * REAL(matrix)[i+j*n] ;
	  REAL(y)[i] = temp ;
	  }
     /*
     ** solve DF'z =b
     */
     for (i=(n-1); i>=0; i--) {
	  if (REAL(matrix)[i+i*n]==0)  REAL(y)[i] =0;
	  else {
	      temp = REAL(y)[i]/REAL(matrix)[i+i*n];
	      for (j= i+1; j<n; j++)
		   temp -= REAL(y)[j]*REAL(matrix)[j+i*n];
	      REAL(y)[i] = temp;
	      }
	  }
	
    return(rank * nonneg);
    }


