;+
;
; Extract a cube slice preserving astro 
;
; input as in HEXTRACT
; waver -> optional wavelength range 
; out*  -> output trimmed variables 
;-


pro ifu_hextract, wave, flux, var, astro, x0, x1, y0, y1, waver=waver,$
                  outwave=outwave, outflux=outflux, outvar=outvar, $
                  outastro=outastro

  ;;trim cube in wave 
  if ~keyword_set(waver) then waver=minmax(wave)

  sz=size(flux)
  wins=where(wave ge waver[0] and wave le waver[1], ninsw)

  mkhdr, tmphdr, flux[*,*,0] 
  putast, tmphdr, astro
  
  for i=0, ninsw-1 do begin

     
     HEXTRACT, flux[*,*,wins[i]], tmphdr, cutflux, cuthead, x0, x1,$
               y0, y1
     
     HEXTRACT, var[*,*,wins[i]], tmphdr, cutvar, tmptmp, x0, x1,$
               y0, y1
     
     ;;init 
     if(i eq 0) then begin
        
        outwave=wave[wins]
        extast, cuthead, outastro
        szcut=size(cutflux)
        outflux=fltarr(szcut[1],szcut[2],ninsw)
        outvar=fltarr(szcut[1],szcut[2],ninsw)
        
     endif
     
     ;;store 
     outflux[*,*,i]=cutflux
     outvar[*,*,i]=cutvar
     
  endfor
  
  

end
