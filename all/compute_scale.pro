;+ 
; NAME:
; compute_scale 
;   Version 1.1
;
; PURPOSE:
;    Apply the scale to images after a aperture have been selected from the first image.
;    Images must be aligned and skysubtracted.
;
; CALLING SEQUENCE:
;   compute_scale,  filename, APPL   
; INPUTS:
;   filename   ASCII file with list of the images to process 
;   APPL       0 to write into file, 1 to apply to images
; RETURNS:
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

pro compute_scale, filename, APPL

       if  N_params() LT 1  then begin 
          print, 'Syntax - ' +$
         'compute_scale,  filename, APPL (v1.1)'
       return
       endif

close, /all


;read the file name

READCOL, filename, name, format="A", /silent

close, /all
 
;select stars
fits=xmrdfits(name[0], 0, header, /fscale, /silent)
skymean=SXPAR(header,'SKYMEAN')
newf=fits+skymean
mwrfits, newf, "tmp.fits", header, /create
markregion, "tmp.fits", REGSTR=reg
spawn, "rm tmp.fits"
;check output
IF(reg.reg_type[1] NE 1 OR reg.nreg EQ 0) THEN BEGIN
print, "No circular regions found.. Skip images", filename[i]
return
ENDIF 

;define array aperture
apercnt=FLTARR(N_ELEMENTS(name),reg.nreg)
scale=FLTARR(N_ELEMENTS(name),reg.nreg)


i=0
WHILE(i LT N_ELEMENTS(name)) DO BEGIN

fits=xmrdfits(name[i], 0, header, /fscale, /silent)

nreg=1
WHILE(nreg LE reg.nreg) DO BEGIN
x0 = reg.regions[nreg,0]  
y0 = reg.regions[nreg,1] 
x1 = reg.regions[nreg,2] 
y1 = reg.regions[nreg,3] 

radius=SQRT((x0-x1)^2+(y1-y0)^2)
index=xpix_circ(x0,y0,radius)

;compute cnt aperture
apercnt[i,nreg-1]=TOTAL(fits[index[0,*],index[1,*]])
nreg=nreg+1
ENDWHILE
i=i+1
ENDWHILE


print, "......"

;determine scale
FINALSCALE=FLTARR(N_ELEMENTS(name))

i=0
WHILE(i LT N_ELEMENTS(name)) DO BEGIN
nreg=1
WHILE(nreg LE reg.nreg) DO BEGIN
scale[i,nreg-1]=apercnt[0,nreg-1]/apercnt[i,nreg-1] 
nreg=nreg+1
ENDWHILE
FINALSCALE[i]=MEDIAN(scale[i,*])
print, "Final scale ",name[i]," ", FINALSCALE[i] 

IF(APPL EQ 1) THEN BEGIN
fits=xmrdfits(name[i], 0, header, /fscale, /silent)
fits=fits*FINALSCALE[i]
commento=STRING("Imaged scaled by ",  FINALSCALE[i])
SXADDPAR, header, "COMMENT", commento
SXADDPAR, header, "IMGSCALE", FINALSCALE[i]
mwrfits, fits,  name[i], header, /create
ENDIF
i=i+1
ENDWHILE

flnam=STRING(filename,".scale")
forprint, name, FINALSCALE, /silent, /NOCOMMENT, TEXTOUT =flnam
 


END
