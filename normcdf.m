function cdf = normcdf (x, mu = 0, sigma = 1)

  if (nargin != 1 && nargin != 3)
    print_usage ();
  endif

  if (! isscalar (mu) || ! isscalar (sigma))
    [retval, x, mu, sigma] = common_size (x, mu, sigma);
    if (retval > 0)
      error ("normcdf: X, MU, and SIGMA must be of common size or scalars");
    endif
  endif

  if (iscomplex (x) || iscomplex (mu) || iscomplex (sigma))
    error ("normcdf: X, MU, and SIGMA must not be complex");
  endif

  if (isa (x, "single") || isa (mu, "single") || isa (sigma, "single"));
    cdf = zeros (size (x), "single");
  else
    cdf = zeros (size (x));
  endif

  if (isscalar (mu) && isscalar (sigma))
    if (isfinite (mu) && (sigma > 0) && (sigma < Inf))
      cdf = stdnormal_cdf ((x - mu) / sigma);
    else
      cdf = NaN (size (x), class (cdf));
    endif
  else
    k = ! isfinite (mu) | !(sigma > 0) | !(sigma < Inf);
    cdf(k) = NaN;

    k = ! k;
    cdf(k) = stdnormal_cdf ((x(k) - mu(k)) ./ sigma(k));
  endif

endfunction