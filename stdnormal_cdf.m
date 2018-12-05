function cdf = stdnormal_cdf (x)

  if (nargin != 1)
    print_usage ();
  endif

  if (iscomplex (x))
    error ("stdnormal_cdf: X must not be complex");
  endif

  cdf = erfc (x / (-sqrt(2))) / 2;

endfunction
