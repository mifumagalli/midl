FUNCTION COADD_REGISTER, z_ref, obj_ref, scifile_shf, shfind $
                         , MIN_SHIFT = MIN_SHIFT

;; Minimum shift is a tenth of a pixel
IF NOT KEYWORD_SET(MIN_SHIFT) THEN MIN_SHIFT = 0.1D
lambda_lya = 1215.67*(1.0+z_ref)
;; Compute the (x,y) position of Lya in the reference image
obj_shf_str = xmrdfits(scifile_shf, 5, /silent)
obj_shf = obj_shf_str[shfind-1L]
ny = n_elements(obj_ref.WAVE_OPT)
;; Compute smooth (S/N)^2 at Lya wavelength
;wave_shf = obj_ref.WAVE_OPT

;dv = 10000.0 ;; med width in km/s
;bkspace = (dv/3.0d5)/alog(10.0d)
;sn2     = (obj_shf.FLUX_OPT^2*obj_shf.IVAR_OPT > 0.0)
;goodpix = where(obj_shf.IVAR_OPT GT 0 AND sn2 GT 0, nzero)
;med_width = round(ny/((max(wave_shf)-min(wave_shf))/bkspace))
;sn2_med = fltarr(ny)
;sn2_med1 = djs_median(sn2[goodpix], width = med_width, boundary = 'reflect')
;sn2_med2 = interpol(sn2_med1, wave_shf[goodpix], wave_shf)
;sig_res = med_width/15L
;nhalf =  long(sig_res)*4L
;xkern = dindgen(2*nhalf+1)-nhalf
;kernel = gauss1(xkern, [0.0, sig_res, 1.0])
;sn2_med = convol(sn2_med2, kernel, /edge_truncate)
;weight = interpol(sn2_med, wave_shf, lambda_lya)  >  0.1D

;print, 'S/N: ', sqrt(weight)
;fwhm = obj_shf.FWHM

;; Find the y-pixel value of Lya
lyapix_y_ref = interpol(findgen(ny), obj_ref.WAVE_OPT, lambda_lya)
lyapix_y_shf = interpol(findgen(ny), obj_shf.WAVE_OPT, lambda_lya)

lyapix_x_ref = interpol(obj_ref.XPOS, findgen(ny), lyapix_y_ref)
lyapix_x_shf = interpol(obj_shf.XPOS, findgen(ny), lyapix_y_shf)

yshift = lyapix_y_ref - lyapix_y_shf 
xshift = lyapix_x_ref - lyapix_x_shf

format = '(%" Measured shift of [%7.2f, %7.2f]")'
splog, xshift, yshift, format = format

IF abs(xshift) LT MIN_SHIFT THEN BEGIN
    xshift = 0.0D
    splog, ' xshift too small, using zero instead'
ENDIF
IF abs(yshift) LT MIN_SHIFT THEN BEGIN
    yshift = 0.0D
    splog, ' yshift too small, using zero instead'
ENDIF

RETURN, [xshift, yshift]
END


