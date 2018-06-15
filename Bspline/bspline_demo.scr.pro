;
; Demonstration script for the initial application
; of radial bspline galaxy-model subtraction.
;
; A. Bolton, 2005 July 06
;

; Supply your own path here:
datadir = '/home/mikifuma/PROGETTI/DLAIMAGING/FINALE/subtract/fitpsf/'
;datadir = ''

; Specify the data filename:
fname = 'DLA1_R_cut.fits'

; Read in the image:
im = mrdfits(datadir+fname)

; Get approximate galaxy positions, which are
; tabulated in an ascii file for now.
posfile = 'quasar_position.txt'
readcol, posfile, xgal, ygal, format='I,I'
ngal = n_elements(xgal)

; To look at the cluster-galaxy positions in the
; image, execute the following:
atv, im
atvplot, xgal, ygal, ps=1

; Set the radius over which to fit, in pixels.
; This is the half-width of the little cutouts that
; we'll extract from the imaage.
rpix = 18

; Get little cutout images:
gstack = fltarr(2*rpix+1,2*rpix+1,ngal)
for i = 0, ngal-1 do gstack[*,*,i] = im[ $
 xgal[i]-rpix:xgal[i]+rpix,ygal[i]-rpix:ygal[i]+rpix]

 
; Don't have error images?  We're probably background
; limited here, so we'll just assume some nominal
; constant variance with which to weight the fit.
invvar = 0.* gstack[*,*,0] + 0.01

; Since the images are drizzled, the RA/Dec grid
; for fitting purposes is trivial.   We'll express
; it in terms of pixels:
ra = reverse(findgen(2*rpix+1)) # replicate(1., 2*rpix+1) - float(rpix)
dec = replicate(1., 2*rpix+1) # findgen(2*rpix+1) - float(rpix)

 
; Perform the centering step.  "rval" controls the radius
; at which the model-based centroiding is performed.
; Note the use of ntheta=[0,-1,1,-2,2], with the -1,1
; entries necessary for allowing for a dipole moment,
; which encodes the centroid shift.
; I set upper=100, lower=100 (rejection thresholds) to
; disable sigma-clipping rejection in the bspline modeling,
; in lieu of actual error values.
ra_cen = fltarr(ngal)
dec_cen = fltarr(ngal)
cen_start = [0., 0.] ; starting guess for center [ra, dec]
rbkpt = 2.* findgen(rpix/2) ; breakpoints every 2 pix in r
ntheta=[0,-1,1,-2,2] ; monopole, dipole, and quadrupole moments.

for i = 0, ngal-1 do begin & $
  newcen = bspline_rtweak(gstack[*,*,i], ra, dec, cen_start, $
   invvar=invvar, rval=2., rbkpt=rbkpt, ntheta=ntheta, ntweak=5, $
   upper=100., lower=100.) & $
  ra_cen[i] = newcen[0] & $
  dec_cen[i] = newcen[1] & $
endfor

; Let's plot up the images and the computed centers, to
; verify that they look plausible:

plotim = gstack[*,*,0]
for i = 1, ngal-1 do plotim = [plotim, gstack[*,*,i]]
atv, plotim
atvplot, rpix + (2*rpix+1)*findgen(ngal) - ra_cen, $
 replicate(rpix, ngal) + dec_cen, ps=1

; Now compute galaxy models assuming these centers.
; You can play around with the degrees of freedom in
; the multipole orders (ntheta=) and the radial
; breakpoint sequence (rbkpt=).  Always include + and -
; ntheta for a given multipole order (those correspond
; to the sin and cos terms).  the orders are:
;  0: monopole
;  -1,1: dipole (exclude if we're centered and symmetric)
;  -2,2: quadrupole (ellipticity)
;  -3,3: trefoil, I guess???  Not sure if this is very useful.
;  -4,4: octupole (cloverleaf)
;  etc. . .
ntheta = [0,-1,1,-2,2]
rbkpt = 2.* findgen(rpix/2)
mstack = 0. * gstack
for i = 0, ngal-1 do begin & $
  delta_ra = ra - ra_cen[i] & $
  delta_dec = dec - dec_cen[i] & $
  r = sqrt(delta_ra^2 + delta_dec^2) & $
  theta = atan(delta_dec, delta_ra) & $
  rset = bspline_radial(r, theta, gstack[*,*,i], invvar=invvar, $
   rbkpt=rbkpt, ntheta=[0,-2,2], yfit=yfit, $
   upper=100, lower=100) & $
  mstack[*,*,i] = yfit & $
endfor

; Finally, have a look at the results:

plotmod = mstack[*,*,0]
for i = 1, ngal-1 do plotmod = [plotmod, mstack[*,*,i]]

atv, [[plotim-plotmod], [plotmod], [plotim]]

; Write this mosaic out to a fits file, for fun:
mwrfits, [[plotim-plotmod], [plotmod], [plotim]], $
 'galmosaic.fits', /create


end
