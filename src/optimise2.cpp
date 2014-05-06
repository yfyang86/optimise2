#include <Rinternals.h>
#include <float.h>		/* for DBL_MAX */
#include <R_ext/Applic.h>	/* for optif9, fdhess */
#include <R_ext/RS.h>	       	/* for Memcpy */
#undef _
#ifdef ENABLE_NLS
#include <libintl.h>
#define _(String) dgettext ("stats", String)
#else
#define _(String) (String)
#endif
#include <R.h>
#include <math.h>
#include <float.h> /* DBL_EPSILON */

#include <Rmath.h>
#include <R_ext/Applic.h>

//#include <RcppGSL.h>
#include <math.h>
//#include <stdlib.h>
#include <string.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_min.h>
#include <gsl/gsl_ieee_utils.h>

#define  GSL_IEEE_MODE double-precision;
struct callinfo {
  SEXP R_fcall;
  SEXP R_env;
} ;

/* fmin(f, xmin, xmax tol) */  
  
extern "C" double fcn1(double x,  void * info1){ 
  struct callinfo * info=(struct callinfo *)info1;
    SEXP s, sx;
    PROTECT(sx = ScalarReal(x));
    SETCADR(info->R_fcall, sx);
    s = eval(info->R_fcall, info->R_env);
    UNPROTECT(1);
    switch(TYPEOF(s)) {
    case INTSXP:
	if (length(s) != 1) goto badvalue;
	if (INTEGER(s)[0] == NA_INTEGER) {
	    warning(_("NA replaced by maximum positive value"));
	    return DBL_MAX;
	}
	else return INTEGER(s)[0];
	break;
    case REALSXP:
	if (length(s) != 1) goto badvalue;
	if (!R_FINITE(REAL(s)[0])) {
	    warning(_("NA/Inf replaced by maximum positive value"));
	    return DBL_MAX;
	}
	else return REAL(s)[0];
	break;
    default:
	goto badvalue;
    }
 badvalue:
    error(_("invalid function value in 'optimize'"));
    return 0;/* for -Wall */
}

extern "C" SEXP optimise2(SEXP call, SEXP op, SEXP args, SEXP rho){
gsl_ieee_env_setup();

try{
   double tol,m,a,b,tracing=0.;
   //m		a	b
   //init	min	max
    SEXP v, res;
    struct callinfo info;

  int status;
  int iter = 0, max_iter = 100;
  const gsl_min_fminimizer_type *T;
  gsl_min_fminimizer *s; 
  gsl_function F;

  /* 
  R part
  */
    args = CDR(args);
    /* the function to be minimized */

    v = CAR(args);
    if (!isFunction(v))
	error(_("attempt to minimize non-function"));
    args = CDR(args);

    /* xmin */

    a = asReal(CAR(args));
    if (!R_FINITE(a))
	error(_("invalid '%s' value"), "xmin");
    args = CDR(args);

    /* xmax */

    b = asReal(CAR(args));
    if (!R_FINITE(b))
	error(_("invalid '%s' value"), "xmax");
    if (a >= b)
	error(_("'xmin' not less than 'xmax'"));
    args = CDR(args);

    /* tol */

    tol = asReal(CAR(args));
    if (!R_FINITE(tol) || tol <= 0.0)
	error(_("invalid '%s' value"), "tol");
    args = CDR(args);
    /*init*/
    m = asReal(CAR(args));
    if (!R_FINITE(m))
	error(_("invalid '%s' value"), "initial");
    args = CDR(args);
    /*tracing*/
    tracing=asReal(CAR(args));
    
    info.R_env = rho;
    PROTECT(info.R_fcall = lang2(v, R_NilValue));
    PROTECT(res = allocVector(REALSXP, 1));
  
  /*gsl*/
  F.function = &fcn1;
  F.params = (void *) &info;
  
     
    if(tracing>.5){
    Rprintf ("%5d [%.7f, %.7f] %.7f %+.7f %.7f\n",
          iter, a, b,
          m, tol , b - a);
    
    }
    
  T =  gsl_min_fminimizer_quad_golden;
  s = gsl_min_fminimizer_alloc (T);
  gsl_set_error_handler_off();
 
    gsl_min_fminimizer_set (s, &F, m, a, b);
 

  if(tracing>.5){
    Rprintf ("using %s method\n",
          gsl_min_fminimizer_name (s));

  Rprintf ("%5s [%9s, %9s] %9s %10s %9s\n",
          "iter", "lower", "upper", "min",
          "err", "err(est)");

  Rprintf ("%5d [%.7f, %.7f] %.7f %+.7f %.7f\n",
          iter, a, b,
          m, m , b - a);
  
  }
  
  
  do
    {
      iter++;
      status = gsl_min_fminimizer_iterate (s);

      m = gsl_min_fminimizer_x_minimum (s);
      a = gsl_min_fminimizer_x_lower (s);
      b = gsl_min_fminimizer_x_upper (s);

      status 
        = gsl_min_test_interval (a, b, tol, 0.000);

      if(tracing>.5){if (status == GSL_SUCCESS)
        Rprintf ("Converged:\n");

      Rprintf ("%5d [%.7f, %.7f] "
              "%.7f %+.7f %.7f\n",
              iter, a, b,
              m, m, b - a);
    }}
  while (status == GSL_CONTINUE && iter < max_iter);
  
    
  gsl_min_fminimizer_free (s);
  REAL(res)[0]=a;
    UNPROTECT(2);
    return res;
  
}catch(...) {
::Rf_error( "c++ exception (unknown reason)" );
}
return R_NilValue; // -Wall
}
