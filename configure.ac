AC_INIT([optimise2], 0.1.0)
# use gsl
AC_PATH_PROG([GSL_CONFIG], [gsl-config])
if test "${GSL_CONFIG}" != ""; then
	GSL_CFLAGS=`${GSL_CONFIG} --cflags`
	GSL_LIBS=`${GSL_CONFIG} --libs-without-cblas`
else
	AC_MSG_ERROR([gsl-config not found, is GSL installed?])
fi
RCPP_LDFLAGS=`Rscript -e 'Rcpp:::LdFlags()'`

AC_SUBST(GSL_CFLAGS)
AC_SUBST(GSL_LIBS)
AC_SUBST(RCPP_LDFLAGS)

AC_OUTPUT(src/Makevars)


