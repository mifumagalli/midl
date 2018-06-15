README.txt

Radial b-spline procedures, such as they are.
Copyright 2005--2007, A. Bolton and S. Burles.

To be regarded as a collection of routines,
not an end-user software package.  Provided
free of charge and without any expressed or
implied warranties.

radial b-splines depend upon code in the idlutils/idlspec2d
IDL tool packages, available through
  http://spectro.princeton.edu/idlspec2d_install.html

Description of the basic idea of radial bsplines can be found in
  Bolton et al. 2006 ApJ, 638, 703
** Please cite this paper in publications that make use of the
** radial bsplines.

Files:

bspline_demo.scr.pro
  (script)
  gives an outline of a deployment of a straightforward
  radial b-spline galaxy modeling application.  Doesn't
  do PSF convolution or anything fancy.  Written to apply
  to multiple galaxy images cut out from a single frame.
  Data not included.

galpos_1252.txt
  goes along w/ above demo script.

bspline_radial.pro
bspline_radial_valu.pro
  actual radial b-spline functions.

bspline_ellip.pro
  advanced b-spline modeling function that can be wrapped
  within MPFIT to do non-linear centering and ellipticity/PA
  solution in the radial/angular coordinates, along with PSF
  convolution and multiple overlapping component modeling.
  Each call to bspline_ellip does a linear b-spline fit to the
  data for the specified non-linear parameter values.

bspline_psf_action.pro
bspline_ellip_valu.pro
  routines related to bspline_ellip

