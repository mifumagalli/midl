;+ 
; NAME:
; uband_cal
;   Version 1.1
;
; PURPOSE:
; For a given Landolt field with SDSS coverage, this procedure compute
; u badn magnitude for the standard stars and rewrite a Landolt list
; in AB mag.
;  
; CALLING SEQUENCE:
;  uband_cal, pgfield    
;
; INPUTS:
;  pgfield  name of the Landolt field to process.
; 
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
 
PRO  uband_cal, pgfield
       if  N_params() LT 1  then begin 
          print, 'Syntax - ' +$
         'uband_cal, pgfield  (v1.1)'
       return
       endif

close, /all

common  colori

fitsname=STRING(pgfield,".fits")
listname=STRING(pgfield,".lst")

fits=mrdfits(fitsname, 0, header, /fscale, /silent)      

READCOL, listname, NAME, RA, DEC, V, BV, UB, VR, RI, VI, n, m, $
         eV, eBV, eUB, eVR, eRI, eVI, FORMAT="A,A,A,F,F,F,F,F,F,F,F,F,F,F,F,F,F",/silent

NAME=STRTRIM(NAME,2)

;read header
MAGZP=SXPAR(header,'MAGZP')    
CTYPE1=SXPAR(header,'CTYPE1')    
CTYPE2=SXPAR(header,'CTYPE2')    
CRVAL1=SXPAR(header,'CRVAL1')    
CRVAL2=SXPAR(header,'CRVAL2')    
CRPIX1=SXPAR(header,'CRPIX1')    
CRPIX2=SXPAR(header,'CRPIX2')    
CDELT1=SXPAR(header,'CDELT1')    
CDELT2=SXPAR(header,'CDELT2')    


;get mag for each star

i=0

Nstar=N_ELEMENTS(RA)

UMAGsdss=FLTARR(Nstar)

WHILE(i LT Nstar) DO BEGIN

RAg=ra_hh2deg(RA[i])
DEg=dec_dd2deg(DEC[i])
 
xcen=(RAg-CRVAL1)/CDELT1+CRPIX1    
ycen=(DEg-CRVAL2)/CDELT2+CRPIX2

fwhm=4.

;set initial guess
IF(i EQ 2) THEN xcen=720
IF(i EQ 2) THEN ycen=1188

;IF(i EQ 0) THEN print, "Initial stars position"
;print, xcen, ycen

CNTRD, fits, xcen, ycen, newx, newy, fwhm, EXTENDBOX=15

;aperture in arcsec
aperture=[10.]
skyape=[12.,16.]
;in pixel
ccdscale=0.3996
aperpix=aperture/ccdscale
skypix=skyape/ccdscale





phpadu=1.
badpix=[-10000,1D8]
X_APER, fits, newx, newy, magscnt, errap, sky, skyerr, $
      phpadu, aperpix, skypix, badpix, /SILENT


;calibrate
UMAGsdss[i]=magscnt+MAGZP-25

i=i+1
ENDWHILE

; compute a bit of color etc...

B=BV+V
U=UB+B
R=V-VR
I=R-RI

;compare my u con U landolt

Uab=U+0.71
Umeanzp=MEAN(UMAGsdss-Uab)
Ustddev=STDDEV(UMAGsdss-Uab)
print, "Mean offset: ", Umeanzp
print, "Dispersion: ", Ustddev

!P.MULTI=[0.,2.,1]

plot, UMAGsdss, Uab, /ynozero, psym=6, ytitle="Landolt (AB)", xtitle="SDSS (AB)" 
oplot, [0,40], [0.,40]

mout=STRING("ZP=",Umeanzp)
xyouts, MIN([UMAGsdss]), MAX([Uab]), mout
sout=STRING("STD=",Ustddev)
xyouts, MIN([UMAGsdss]), MAX([Uab])-0.2, sout


;derive ab mag and colors

u_AB=UMAGsdss
B_AB=B-0.15
V_AB=V
R_AB=R+0.199
I_AB=I+0.454

BVab=B_AB-V_AB
uBab=u_AB-B_AB
VRab=V_AB-R_AB
RIab=R_AB-I_AB
VIab=V_AB-I_AB

;add rms to color
neweUB=SQRT(eUB^2+Ustddev^2)


plot, BVab, BV, /ynozero, psym=6, xtitle="(AB) colors", ytitle="Landolt colors",yrange=[-2,2],xrange=[-2,2]
oplot, uBab, UB,  psym=6, color=rosso
oplot, VRab, VR,  psym=6, color=blu
oplot, RIab, RI,  psym=6, color=giallo
oplot, VIab, VI,  psym=6, color=verde
oplot, [-5,5], [-5,5] 

!P.MULTI=0  


;print plot

plotname=STRING(pgfield,"_cal.ps")
SET_PLOT, "PS"
DEVICE, filename=plotname, /color
!P.MULTI=[0.,2.,1]
plot, UMAGsdss, Uab, /ynozero, psym=6, ytitle="Landolt (AB)", xtitle="SDSS (AB)" 
oplot, [0,40], [0.,40]
xyouts, MIN([UMAGsdss]), MAX([Uab]), mout
xyouts, MIN([UMAGsdss]), MAX([Uab])-0.2, sout
plot, BVab, BV, /ynozero, psym=6, xtitle="(AB) colors", ytitle="Landolt colors",yrange=[-2,2],xrange=[-2,2]
oplot, uBab, UB,  psym=6, color=rosso
oplot, VRab, VR,  psym=6, color=blu
oplot, RIab, RI,  psym=6, color=giallo
oplot, VIab, VI,  psym=6, color=verde
oplot, [-5,5], [-5,5] 
!P.MULTI=0  
DEVICE, /CLOSE
SET_PLOT, "X"



; rewrite landolt file

newname=STRING(pgfield,"AB.lst")
nomicoord=STRARR(N_ELEMENTS(NAME))

FOR i=0, N_ELEMENTS(NAME)-1 DO nomicoord[i]=STRING(NAME[i]," ",RA[i]," ",DEC[i])





FORPRINT, nomicoord,V_AB, BVab, uBab, VRab, RIab, VIab, n, m, $
         eV, eBV, neweUB, eVR, eRI, eVI, FORMAT="A,F,F,F,F,F,F,F,F,F,F,F,F,F,F",$
	 /silent, /NoCOMMENT, TEXTOUT=newname 


;update structure

stru=mrdfits("/home/mikifuma/idl/xidl/IMG/Photometry/Lists/nlandoltAB.fits",1)

stru.NAME=STRTRIM(stru.NAME,2)

i=0
WHILE(i LT Nstar) DO BEGIN
index=WHERE(stru.NAME EQ NAME[i])
stru[index].BV=BVab[i]
stru[index].UB=uBab[i]
stru[index].VR=VRab[i]
stru[index].RI=RIab[i]
stru[index].VI=VIab[i]
stru[index].SIG_UB=neweUB[i]
i=i+1
ENDWHILE


mwrfits, stru, "/home/mikifuma/idl/xidl/IMG/Photometry/Lists/nlandoltAB.fits", /CREATE, /SILENT


end


