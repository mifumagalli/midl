;+ 
; NAME:
; final_trim
;   Version 1.1
;
; PURPOSE:
; Normalise images in count/sec using exp field in X structure. Update header and structure. 
;
; CALLING SEQUENCE:
;   final_trim, imgin, imgout, ROT=, X0=, Y0=, X1=, Y1=   
;
; INPUTS:
;  imgin        image to trim
;  imgout       name for rinal images
;  ROT          0: none, 1: +90 counterclock
;  X0,Y0,X1,Y1  edge of final images
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
  
PRO  final_trim, imgin, imgout, ROT=rot, X0=x0, Y0=y0, X1=x1, Y1=y1
  
       if  N_params() LT 1  then begin 
          print, 'Syntax - ' +$
         'final_trim, imgin, imgout, ROT=, X0=, Y0=, X1=, Y1= (v1.1)'
       return
       endif

close, /all


fits=xmrdfits(imgin, 0, header, /fscale, /silent)

fits=fits[x0:x1,y0:y1]

IF(ROT EQ 1) THEN fits=rotate(fits,1)

mwrfits, fits,  imgout, header, /create

end
