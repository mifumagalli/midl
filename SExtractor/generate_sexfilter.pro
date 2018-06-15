;+
;
; Generate a filter file for sextractor
; Currently supported only gaussian filter
;
;
;
; fwhm   the FWHM in pixel of the gaussian
; size   the size of the box
;
;
;-


pro generate_sexfilter, fwhm, size

;;make int
size=fix(size)

;;make size odd
if(size mod 2 eq 0) then size=size+1

;;make the gaussian
gauss2d,size,size,floor(0.5*size),floor(0.5*size),fwhm,gaussian

;;open unit
openw, 1, 'gauss_'+string(fwhm,format='(F3.1)')+'_'+rstring(size)+'x'+rstring(size)+'.conv'


stop

;;print the file
printf, 1, 'CONV NORM'
printf, 1, '#'+rstring(size)+'x'+rstring(size)+$
  ' convolution mask of a gaussian PSF with FWHM = '+string(fwhm,format='(F3.1)')+' pixels.'
printf, 1, gaussian


;;print to screen
print, '-------------------'
print, gaussian
print, '-------------------'


close, /all





end
