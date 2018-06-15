;+
; NAME:     KUBEVIZ_ZSCALE_RANGE
;
; PURPOSE: Computes the best levels cuts for displaying astronomical
;	   images using IRAF zscale algorithm.
;
;
; DESCRIPTION: Returns the best low_cut and high cut values for the
;		visualization of an image, on the base of the IRAF
;		Zscale algorithm.
;
; CATEGORY: Image Visualization.
;
; CALLING SEQUENCE: Result = ZSCALE_RANGE(image, contrast, POINTS=POINTS)
;
; INPUTS: Image:    1-D or 2-D array.
;         Contrast: Constrast for the Zscale algorithm; default = 1.
;
; KEYWORDS:
;
;       POINTS: Number of sample points to use for the determination
;       of  the range; default = 100 points.
;
; OUTPUTS: 2-element vector of the form [low_cut, high-cut].
;
; PROCEDURE:  Makes a linear fit to the sorted values of the sample
;             points and takes the minumum and maximum values of the fit.
;
; MODIFICATION HISTORY:
;
;       Feb 2004 - Gianluca Li Causi, INAF - Rome Astronomical Observatory
;				licausi@mporzio.astro.it
;				http://www.mporzio.astro.it/~licausi/
;       Oct 2015 - Matteo Fossati renamed into kubeviz_zscale_range to avoid
;                                 conflicts
;
;-

FUNCTION kubeviz_zscale_range, im, contrast, points=points

IF n_elements(contrast) EQ 0 THEN contrast = 1.

IF contrast LE 0 THEN contrast = 1.

IF min(im, /nan) EQ max(im, /nan) THEN RETURN, [min(im, /nan), max(im, /nan)]

IF n_elements(points) EQ 0 THEN points = 100

siz = size(im)


INIZIO:
CASE siz[0] OF
    
    1: BEGIN
        lx = siz[1]
        np = points < lx	;numero di punti per la griglia
        xx = findgen(np) * (lx-1)/(np-1)
        imp = im[xx]
    END
    
    2: BEGIN
        lx = siz[1]
        ly = siz[2]
        
        np = points < min([lx, ly]) ;numero di punti per la griglia
        	
        xx = findgen(np) # replicate(1, np) * (lx-1)/(np-1)
        yy =  replicate(1, np) # findgen(np) * (ly-1)/(np-1)
        
        imp = im[xx,yy]

    END
    
ELSE:
    
ENDCASE

IF min(imp, /nan) EQ max(imp, /nan) THEN BEGIN
	points = points * 5.
	GOTO, INIZIO
ENDIF

s = sort(imp)
intens = imp[s]

;Remove nans
valid = where(finite(intens) eq 1, Nfinite)
intens = intens[valid]


;x = findgen(np^siz[0])
x = findgen(Nfinite)

;plot, x, intens

; linear fit
niter = 3
soglia = 3.
coef = linfit(x, intens, yfit=yfit)
FOR i = 1, niter-1 DO BEGIN
	diff = ABS(intens - yfit)
	IF n_elements(diff) LE 1 THEN BREAK
	sig = stdev(diff)
	ok = where(diff LT soglia*sig, count)
	IF count LE 1 THEN BREAK
	x = x[ok]
	intens = intens[ok]
	coef = linfit(x, intens, yfit=yfit)
ENDFOR

zmed = (min(yfit) + max(yfit))/2.
xmed = (min(x) + max(x))/2.
zmax = zmed + coef[1] * (max(x)-xmed) / contrast < max(im, /nan)
zmin = zmed + coef[1] * (min(x)-xmed) / contrast > min(im, /nan)

IF zmin EQ zmax THEN BEGIN
	zmin = min(im, /nan)
	zmax = max(im, /nan)
ENDIF

RETURN, [zmin, zmax]

END
