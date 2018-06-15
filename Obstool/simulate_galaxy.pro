;+
; NAME:
;       SIMULATE_GALAXY
; PURPOSE:
;       Create an array within IDL with an oversampled image of a
;       sersic profile with specified galaxy parameters (centered
;       apart from a slight offset of subpixels). It shows the
;       BASIC step of the extensive simulations used for the
;       code-testing of GALFIT and GIM2D shown in the paper Haeussler
;       et al. 2006
; EXPLANATION:
;       This procedure takes the imput parameters, simulates a highly
;       oversampled galaxy and returns this profile as an image.
;
; CALLING SEQUENCE:
;       simulate_galaxy,image, image, mag, re, AR, PA, n, offx, offy
; EXAMPLE
;  simulate_galaxy, image, 21., 15., 0.5, 90, 1., 0.0, 0.0, 24.862, rawimage='/data/images/s9z99A_sci_flx.fits'
;
; INPUTS:
;       mag - magnitude of the simulated galaxy
;       radius - size of the galaxy in pixels (the code was tested
;                from 1 pixel to 500 pixels)
;       ar - axis ratio of the galaxy. Be careful not to use 0, that
;            might cause severe problems
;       pa - position angle (counterclockwise, 0 is horizontal!)
;       sersic - sersic index of the profile
;       offsetx, offsety - offset within pixel in X and Y direction
;               these numbers are tricky. Currently numbers with
;               precision 0.1 are possible (any number in steps of
;               0.1). Due to different definitions where in the pixels
;               (0,0) is defined, these number varies for every program.
;               Higher precsion offsets are in principle possible, but
;               not needed for our purpose.
;
;               This shows two examples how/where (0,0) can be defined
;               within a pixel
;
;                 e.g. GALFIT            e.g.THIS CODE
;                _____________          _____________
;               |             |        |             |
;               |             |        |             |
;               |             |        |             |
;               |             |      0 |             |
;               |             |        |             |
;               |             |        |             |
;             0 |_____________|        |_____________|
;               0                             0
;      
;               E.g. if you set a galaxy on offsetx/y (0.0,0.0), GALFIT
;               will tell you it is centered at (0.5, 0.5)
;       zeropoint - 1 cnt/sec ~ magnitudes (e.g. 24.862 represents z-Band zeropoint of the HST-ACS camera)
;       exptime - you can specify exposure time here or give a raw
;                 image where the exposure time should be read from the header 
;       rawimage - the image the profile should be set into in the end. Used only to read
;                  in the image header for zeropoint calibration (Keyword EXPTIME)
;
; OUTPUTS:
;       image - image containing the galaxy (SERSIC) profile
;
; NOTES:
;       currently this procedure only creates SERSIC profiles. Changing
;       it to create any other profile should be easy and is left to
;       the user
;       This routine returns a noisefree, PSF-free profile. adding
;       this image to a bigger frame (a loop for several profeils) is
;       trivial and left to the user.
;       Also, adding noise and convolving the image with a PSF (and
;       maybe adding it to a real backgroud image) is trivial and not
;       shown here
;
;       The profiles are oversampled in the very inner regions by a
;       factor of 100, middle areas by factor of 10. The image is not
;       completely oversampled because of memory reasons. The needed
;       array would be too big to be handled by IDL on normal computers
;       That is why oversampling areas have a maximum size.
;
; REVISION HISTORY:
;       Written  B. Haeussler upto August 2006
;-

PRO get_kappa, n, k
  COMMON sersickappa, sersickapp
  sersickapp = double(n)
  k = zbrent(0., 20., Func_name = 'kappa_funct',MAX_ITERATIONS = 50)
END

FUNCTION kappa_funct, kappa
  COMMON sersickappa, sersickapp
  sersn = double(sersickapp)
  f = GAMMA(2*sersn)-2*IGAMMA(2*sersn, kappa)*Gamma(2*sersn)
  return, f
END



PRO simulate_galaxy, image, mag, radius, ar, pa, sersic, offsetx, offsety, $
                     zeropoint, exptime = exptime,  rawimage = rawimage

; change number formats to double for higher precision calculations
  mag = double(mag)
  ar = double(ar)
  radius = double(radius)
  pa = double(pa)
  sersic = double(sersic)

; get kappa for the sersic profile
  get_kappa, sersic, kappa
  kappa = double(kappa)

; get correct zeropoint corrected for exposure time
  IF NOT keyword_set(exptime) THEN BEGIN
    IF NOT keyword_set(rawimage) THEN BEGIN
      print, 'You have to specify either exposure time or raw image to read header keyword "EXPTIME"'
      stop
    ENDIF
    head = HEADFITS(rawimage, /silent)
    exptime = sxpar(head, 'EXPTIME')
  ENDIF 
  IF NOT keyword_set(exptime) THEN BEGIN
    print, 'no exposure time specified or found in image header' 
    stop
  ENDIF 
  zeropoint = zeropoint+2.5*alog10(exptime)

; transfer magnitude to total flux (counts)
  flux = 10^((zeropoint - mag)/2.5)

; estimate inner part first, oversample by 100
; set size of inner postage stamp (min for precision, max for memory)
  sizeinn = long(200*round(radius/float(ar))  < 2000 > 1000)   ;max 20x20 speed and memory, min output 10x10
; scale up parameters
  radiusinn = radius*100.*ar
  offsetxinn = offsetx*100.
  offsetyinn = offsety*100.
  simposxinn = sizeinn/2+offsetxinn-0.5
  simposyinn = sizeinn/2+offsetyinn-0.5
; simulate profile
  dist_ellipse,iminn,[sizeinn, sizeinn],simposxinn, simposyinn,ar,pa, /double
  gamma = Gamma(2*sersic)
  iminn = flux*(ar/(2*!pi*radiusinn^2*exp(kappa)*sersic*kappa^(-2*sersic)*gamma))*$
    exp(-kappa*((iminn/radiusinn)^(1/sersic)-1))
; rebin image to have correct size
  iminn = frebin(iminn, sizeinn/100., sizeinn/100., /total)

; estimate middle part then, oversample by 10
; set size of middle postage stamp (min for precision, max for memory)
  sizemid = long(20*round(radius/float(ar))  < 2000 > 1000)   ;max 200x200 speed and memory, min output 100x100
; scale up parameters
  radiusmid = radius*10.*ar
  offsetxmid = offsetx*10.
  offsetymid = offsety*10.
  simposxmid = sizemid/2+offsetxmid-0.5
  simposymid = sizemid/2+offsetymid-0.5
; simulate profile
  dist_ellipse,immid,[sizemid, sizemid],simposxmid, simposymid,ar,pa, /double
  gamma = Gamma(2*sersic)
  immid = flux*(ar/(2*!pi*radiusmid^2*exp(kappa)*sersic*kappa^(-2*sersic)*gamma))*$
    exp(-kappa*((immid/radiusmid)^(1/sersic)-1))
; rebin image to have correct size
  immid = frebin(immid, sizemid/10., sizemid/10., /total)

; simulate outer region
  size = long(100*round(radius/float(ar)) < 4000 > 500)    ; max 4000x4000 speed and memory, min output 500 x 500
;  size = 500
  radius = radius*ar
  simposx = size/2+offsetx-0.5
  simposy = size/2+offsety-0.5
; simulate profile
  dist_ellipse,im,[size, size],simposx, simposy,ar,pa, /double
  gamma = Gamma(2*sersic)
  im = flux*(ar/(2*!pi*radius^2*exp(kappa)*sersic*kappa^(-2*sersic)*gamma))*$
    exp(-kappa*((im/radius)^(1/sersic)-1))


; create final image and put in the three different areas
  image = temporary(im)

; put in middle part
  image[n_elements(image[*,1])/2.-n_elements(immid[*,1])/2.: $
        n_elements(image[*,1])/2.+n_elements(immid[*,1])/2.-1, $
        n_elements(image[1,*])/2.-n_elements(immid[*,1])/2.: $
        n_elements(image[1,*])/2.+n_elements(immid[*,1])/2.-1] = immid

; put in inner part
  image[n_elements(image[*,1])/2.-n_elements(iminn[*,1])/2.: $
        n_elements(image[*,1])/2.+n_elements(iminn[*,1])/2.-1, $
        n_elements(image[1,*])/2.-n_elements(iminn[*,1])/2.: $
        n_elements(image[1,*])/2.+n_elements(iminn[*,1])/2.-1] = iminn

  undefine, immid, iminn, im

end
