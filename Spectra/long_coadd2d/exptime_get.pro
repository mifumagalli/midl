;+
;PURPOSE
;	to get the exposure time out of a header for many telescopes
;	with many types of header layours
;SYNTAX
;	exptime=exptime_get(hdr)
;INPUTS
;	hdr: the header of the fits file
;OUTPUTS
;	exptime: the exposure time in seconds
;NOTES
;	pulled from long_extinct.pro
;
;Written by R. da Silva, 7-27-09, UCSC
;-

FUNCTION exptime_get, scihdr

  compile_opt strictarr

longslit_dir = getenv('LONGSLIT_DIR')
; read in atmospheric extinction data
IF strmatch(sxpar(scihdr, 'TELESCOP'), 'Keck*') THEN BEGIN
    extinctfile = longslit_dir + '/calib/extinction/mkoextinct.dat'
    exptime = double(sxpar(scihdr, 'ELAPTIME')) 
    airmass = double(sxpar(scihdr, 'AIRMASS'))
ENDIF ELSE IF strmatch(sxpar(scihdr, 'OBSERVAT'), 'Gemini*') THEN BEGIN
    if strmatch(sxpar(scihdr, 'TELESCOP'), 'Gemini-North') THEN BEGIN
        extinctfile = longslit_dir + '/calib/extinction/mkoextinct.dat'
    ENDIF ELSE IF strmatch(sxpar(scihdr, 'TELESCOP'), 'Gemini-South') THEN $
        extinctfile = longslit_dir + '/calib/extinction/ctioextinct.dat'
    exptime = double(sxpar(scihdr, 'ELAPSED')) > $ 
              double(sxpar(scihdr, 'EXPTIME')) 
    airmass = double(sxpar(scihdr, 'AIRMASS'))
ENDIF ELSE IF strmatch(sxpar(scihdr, 'DETECTOR'), 'mars*') OR $
   strmatch(sxpar(scihdr, 'DETECTOR'), '*ccd34*') OR $ 
   strmatch(sxpar(scihdr, 'TELESCOP'), 'mmt*') OR  $
   strmatch(sxpar(scihdr, 'TELESCOP'), '*kp4m*') OR $
   strmatch(sxpar(scihdr, 'TELESCOP'), '*3.5m*') $
THEN BEGIN
    extinctfile = longslit_dir + '/calib/extinction/kpnoextinct.dat'
    exptime = double(sxpar(scihdr, 'OBSTIME')) > double(sxpar(scihdr, 'EXPTIME'))
    airmass = double(sxpar(scihdr, 'AIRMASS'))
ENDIF ELSE if (stregex(sxpar(scihdr,'INSTRUME'),'.*kast.*',$
                       /boolean,/fold_case) eq 1) or $
  (stregex(sxpar(scihdr,'VERSION'),'kast*', /bool, /fold_case) EQ 1) then begin
    extinctfile = longslit_dir + '/calib/extinction/mthamextinct.dat'
    exptime = double(sxpar(scihdr, 'EXPOSURE') > sxpar(scihdr,'EXPTIME'))
    ;; Calculate airmass
    ras = sxpar(scihdr, 'RA')
    decs = sxpar(scihdr, 'DEC')
    x_radec, ras, decs, rad, decd
    hangl = float(strsplit(sxpar(scihdr,'HA'),':',/extrac))
    if hangl[0] LT 0. then hangl = hangl[0]-hangl[1]/60. $
    else hangl = hangl[0]+hangl[1]/60. 
    airmass = airmass(40., decd, hangl)
ENDIF ELSE IF  strmatch(strtrim(sxpar(scihdr, 'TELID'), 2), '200') THEN BEGIN
    print, 'Warning:  Using KPNO for Palomar!'
    extinctfile = longslit_dir + '/calib/extinction/kpnoextinct.dat'
    exptime = double(sxpar(scihdr, 'EXPTIME'))
    airmass = double(sxpar(scihdr, 'AIRMASS'))
ENDIF ELSE message, 'Telescope not found'
IF NOT KEYWORD_SET(exptime) and not keyword_set(NOTIME) THEN $
  message, 'Could not find exposure time'

return, exptime

END
