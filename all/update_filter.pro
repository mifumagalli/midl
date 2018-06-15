;+ 
; NAME:
; update_filter
;   Version 1.1
;
; PURPOSE:
; Update the U filter for images in the X structure for imaging reduction.
;
; CALLING SEQUENCE:
; update_filter, struct   
;
; INPUTS:
;  stuct X structure for imagin reduction.
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
 
PRO update_filter, struct

       if  N_params() LT 1  then begin 
          print, 'Syntax - ' +$
         'update_filter, struct (v1.1)'
       return
       endif

close, /all


blue = where(STRPOS(struct.img_root,"blue") NE -1)
struct[blue].filter='U'	
print, "New filters ", struct.filter

; Resave the updated structure so that you can pick up where you left off easily.
mwrfits, struct, 'struct.fits', /create
END
