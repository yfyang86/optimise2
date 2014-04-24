optimise2
=========

The optimise(stat) in R sometimes fails to converge to global optimization point. In fact, in GSL, most of these cases will report a warining message. The rule in this packages is gsl_min_fminimizer_quad_golden with more than one guess.

    Dependencies: Rcpp RcppGSL

In fact, this package may only depend on libgsl. I fail to configure it on windows if Cygwin is not used.


