;+ 
; NAME:
; shift_image
;   Version 1.1
;
; PURPOSE:
; Shift images for a known offset
;
; CALLING SEQUENCE:
;   shift_image, img, offset    
;
; INPUTS:
;  lista list of images
;  img      file with a list of images to shift  
;  offset   file with a list of offset x and y
; 
; RETURNS:
;
; OUTPUTS:
;  
;
; OPTIONAL OUTPUTS:
;
; COMMENTS:
;
; EXAMPLES:
;
; PROCEDURES/FUNCTIONS CALLED:
;  
;  
;
; REVISION HISTORY:
;   11-Nov-2008 Written by MF
;-
;------------------------------------------------------------------------------


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
 
PRO shift_image, img, offset    

       if  N_params() LT 1  then begin 
          print, 'Syntax - ' +$
         'shift_image, img, offset    (v1.1)'
       return
       endif

close, /all



READCOL, img, filename, format="A", /silent
READCOL, offset, xshif, yshif, format="F,F", /silent

i=0
WHILE(i LT N_ELEMENTS(filename)) DO BEGIN

fits=xmrdfits(filename[i], 0, header, /fscale, /silent)

print, "Shifting  ", filename[i]
shfits=SHIFT(fits,xshif[i],yshif[i])

;write file

name=STRING("s_",filename[i])
mwrfits, shfits, name, header, /create

i=i+1
ENDWHILE


END
