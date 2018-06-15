;+ 
; NAME:
; create_mask 
;   Version 1.1
;
; PURPOSE:
; Create a mask by dividing two images for each image to be combined. 
; Read IMGSCALE (apply_scale) SKYMEAN and SKYSIGMA (marksky and skyfit). Update the header.
;
; CALLING SEQUENCE:
;   create_mask, imgtomask, other,  N    
;
; INPUTS:
;  imgtomask image for which you want a mask
;  other     other image to build the mask
;  N         sigma value to reject.  
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
 
PRO create_mask, imgtomask, otherimg, N 

       if  N_params() LT 1  then begin 
          print, 'Syntax - ' +$
         'create_mask, imgtomask, otherimg, N  (v1.1)'
       return
       endif

close, /all





fits1=xmrdfits(imgtomask, 0, header, /fscale, /silent)
fits2=xmrdfits(otherimg, 0, head, /fscale, /silent)



 
SCL1=SXPAR(Header,'IMGSCALE')
SCL2=SXPAR(Head,'IMGSCALE')
BCK1=SXPAR(Header,'SKYMEAN')
BCK2=SXPAR(Head,'SKYMEAN')
SKYS1=SXPAR(Header,'SKYSIGMA')
SKYS2=SXPAR(Head,'SKYSIGMA')

IF(SCL1  LE 0.) THEN SCL1=1.
IF(SCL2  LE 0.) THEN SCL2=1.
IF(BCK1  LE 0.) THEN BCK1=1.
IF(BCK2  LE 0.) THEN BCK2=1.
IF(SKYS1 LE 0.) THEN SKYS1=1.
IF(SKYS2 LE 0.) THEN SKYS2=1.





;reset the level background
ratio=fits1
nonzero=where(fits2 NE 0.)
ratio[nonzero]=(fits1[nonzero]/SCL1+BCK1)/(fits2[nonzero]/SCL2+BCK2)
scale1=MEAN(fits1[nonzero]/SCL1+BCK1)
scale2=MEAN(fits2[nonzero]/SCL2+BCK2)
ratio[nonzero]=ratio[nonzero]*(scale2/scale1)
 
 
 
;write mask
mask=ratio-ratio
soglia=1+N*SQRT(SKYS1^2+SKYS2^2)
cosmic=where(ratio GT soglia, COMPLEMENT=data, numb)
IF(numb GT 0.) THEN mask[cosmic]=1.
name=STRING("mask_",imgtomask)

;mask[*,*]=0
;mask[806,*]=1.
;mask[810,*]=1.


mwrfits, mask, name, hh, /create

print, "Created mask ", name

;update the images header with the BPM
SXADDPAR, header, "BPM", name

mwrfits, fits1,  imgtomask, header, /create

END
