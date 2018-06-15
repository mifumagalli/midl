;+ 
; NAME:
; combweigh, lista
;   Version 1.1
;
; PURPOSE:
;  Compute and add to the header a weigh. It reads from the header IMGSCALE (apply_scale.pro) 
;  and SKYSIGMA (skystat). Produce a file listain.wgh for imcombine.
;
; CALLING SEQUENCE:
;  combweigh, lista  
; INPUTS:
;   lista      ASCII file name of images to process 
;
; RETURNS:
;
;
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

pro  combweigh, lista

       if  N_params() LT 1  then begin 
          print, 'Syntax - ' +$
         'combweigh, lista (v1.1)'
       return
       endif

close, /all



;read the file name

READCOL, lista, name, format="A", /silent

PESO=FLTARR(N_ELEMENTS(name))
i=0
WHILE(i LT N_ELEMENTS(name)) DO BEGIN
fits=xmrdfits(name[i], 0, header, /fscale, /silent)
SIGMA=SXPAR(Header,'SKYSIGMA')
PESO[i]=1./(SIGMA)^2
SXADDPAR, header, "COMBWEIGH", PESO[i]
mwrfits, fits,  name[i], header, /create
print, name[i], PESO[i]
i=i+1
ENDWHILE
newfile=STRING(lista,".wgh")
forprint, PESO, /SILENT, TEXTOUT =newfile,/NoCOMMENT
print, "Created file: ", newfile 

END
