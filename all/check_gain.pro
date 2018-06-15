    
PRO  check_gain, img, chip
  
       if  N_params() LT 1  then begin 
          print, 'Syntax - ' +$
         'check_gain, img, chip (v1.1)'
       return
       endif

close, /all
  
   

; set gain for blue images
 
 
IF (chip EQ 1) THEN  BEGIN
fits=xmrdfits(img, 0, header, /fscale, /silent)
fits[3277:4619,*]=1.66*fits[3277:4619,*]
fits[2252:3276,*]=1.66*fits[2252:3276,*]
fits[0:2251,*]=1.55*fits[0:2251,*]
;notte2 1.68, 1.64, 1.55
name=STRING("testg_",img)
mwrfits, fits,  name, header, /create
print, "Gain applied to ", name
xatv, name
ENDIF

; set gain for red images

IF (chip EQ 2) THEN  BEGIN
fits=xmrdfits(img, 0, header, /fscale, /silent)
fits[0:1063,*]=1.99*fits[0:1063,*] 
fits[1064:2247,*]=2.15*fits[1064:2247,*] 
;notte2 2.00, 2.13
name=STRING("testg_",img)
mwrfits, fits,  name, header, /create
print, "Gain applied to ", name
xatv, name
ENDIF

end
