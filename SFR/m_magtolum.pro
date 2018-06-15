;  Convert AB mag to luminosity

;  inmag -- Magnitude
;  H0 -- Hubble constant
;  FLUX --> if selected, return flux
;  zred --> redshift of the object. If not set assume Distance=10pc
;  kcor --> if set, simple k-correction for flat sed performed

function m_magtolum, inmag, h0=h0, FLUX=flux, ZRED=ZRED, KCOR=KCOR

  if  N_params() LT 1 then begin 
      print,'Syntax - ' + $
        'lum = x_magtolum(absmag, [h0]) (v1.0)'
      return, -1
  endif 

  ;; Convert mag according to h0
  if keyword_set(h0) then absmag = inmag - 5*alog10(h0) $
  else absmag = inmag


  ;; Flux
  flu = 10^(-1.d*(absmag + 48.6)/2.5)

  ;; Lum
  c = x_constants()
  
  IF keyword_set(ZRED) THEN BEGIN
  distance_calculator, ZRED, Dist, /lum, /DEF
  lum = flu * 4 * !pi * (Dist*1D6*c.pc)^2
  ENDIF 
  
  IF keyword_set(KCOR) THEN lum = lum/(1+ZRED)
  
  IF ~keyword_set(ZRED) THEN lum = flu * 4 * !pi * (10*c.pc)^2
 
  IF keyword_set(FLUX) THEN  return, Flu  ELSE  return, lum 
 
end


