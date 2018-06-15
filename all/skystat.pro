;+ 
; NAME:
; skystat
;   Version 1.1
;
; PURPOSE:
; Compute the sky MEAN and STD.  If  SUBT=1, subtract.
;
; CALLING SEQUENCE:
;   skystat, img, SUBT     
;
; INPUTS:
;  lista list of images
;  img      file with a list of images to pocess  
;  SUBT     0 to evaluate only and 1 to subtract
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
 
PRO skystat, img, SUBT      

       if  N_params() LT 1  then begin 
          print, 'Syntax - ' +$
         'skystat, img, SUBT      (v1.1)'
       return
       endif

close, /all



READCOL, img, filename, format="A", /silent
 
i=0
WHILE(i LT N_ELEMENTS(filename)) DO BEGIN

fits=xmrdfits(filename[i], 0, header, /fscale, /silent)

;select sky region
markregion, filename[i], REGSTR=reg

;evaluate mean, sigma

IF(reg.reg_type[0] NE 0 OR reg.nreg EQ 0) THEN BEGIN
print, "No square regions found.. Skip image", filename[i]
ENDIF 

IF(reg.reg_type[0] EQ 0 AND reg.nreg NE 0) THEN BEGIN

meansky=FLTARR(reg.nreg)
sigmasky=FLTARR(reg.nreg)
edge=STRARR(reg.nreg)

nreg=1
WHILE(nreg LE reg.nreg) DO BEGIN
x0 = reg.regions[nreg,0]  
y0 = reg.regions[nreg,1] 
x1 = reg.regions[nreg,2] 
y1 = reg.regions[nreg,3] 
     

;compute first mean and stddev
meansky[nreg-1]=MEAN(fits[x0:x1,y0:y1])
sigmasky[nreg-1]=STDDEV(fits[x0:x1,y0:y1]) 

iter=0
WHILE(iter LT 5) DO BEGIN
;iterate to reject bright objects
newreg=fits[x0:x1,y0:y1]
index=where(newreg LT (meansky[nreg-1]+3*sigmasky[nreg-1]),cnt)
;compute new sigma after rejection
sigmasky[nreg-1]=STDDEV(newreg[index]) 
;compute new mean after rejection
meansky[nreg-1]=MEAN(newreg[index])
iter=iter+1
ENDWHILE


edge[nreg-1]=STRING(x0,",",x1,",",y0,",",y1)


IF(nreg EQ 1) THEN print, "................ Image ", filename[i]
print, "Region ", nreg 
print, "Edges ", edge[nreg-1]
print, "Mean ", meansky[nreg-1] 
print, "Sigma ", sigmasky[nreg-1] 

s=STRING(nreg)
key=STRING("SKYBOX",STRTRIM(s,2))
SXADDPAR, header, key, edge[nreg-1]
 
nreg=nreg+1
ENDWHILE

allmean=MEAN(meansky)
allsigma=MEAN(sigmasky)

print, "Overall mean sky ", allmean 
print, "Overall std sky ", allsigma

SXADDPAR, header, "SKYMEAN", allmean 
SXADDPAR, header, "SKYSIGMA", allsigma

IF(SUBT EQ 1) THEN fits=fits-allmean 

mwrfits, fits,  filename[i], header, /create

ENDIF


i=i+1
ENDWHILE



END
