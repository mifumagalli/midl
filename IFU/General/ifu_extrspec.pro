;+
; 
; Take a data cube and extract a spectrum in a given aperture 
; wave       input wave 
; cube       input cube  
; variance   associated variance 
; center     array with center (ra,dec) 
; aperture   dimatere of parture in "
; outspc     output spectrum in units of cube 
; outvar     output variance
; astro      the astro structure 
; write      name of the output name 
; show       display spectrum 
;-

pro ifu_extrspec, wave, cube, var, center=center, aperture=aperture, $
                  outspc=outspc, outvar=outvar, astro=astro, write=write, show=show
  
  
  ;;find corresponding xy location 
  x_radec, center[0],  center[1], rag, deg
  ad2xy, rag, deg, astro, x0, y0
  
  ;;size in arcsec 
  psca=astro.cd[1,1]*3600.
  radius=aperture/psca/2.
  
  ;;find enclosed pixels 
  sz=size(cube)
  pix=mpix_circ(x0, y0, radius, XSIZE=sz[1], YSIZE=sz[2])  

  outspc=fltarr(sz[3])
  outvar=fltarr(sz[3])
  outwgt=fltarr(sz[3])
  
  for cc=0, sz[3]-1 do begin
 
     slice=cube[*,*,cc]
     slice_var=var[*,*,cc]


     ;;find nans and kill them in variance 
     nanfl=where(slice[pix] ne slice[pix], nnnan)
     mask=fltarr(sz[1],sz[2])+1.
     if(nnnan gt 0) then  mask[nanfl]=0.
   
     ;;add masking nans 
     outspc[cc]=total(mk_finite(slice[pix]/slice_var[pix])*mask)
     outwgt[cc]=total(mk_finite(1./slice_var[pix])*mask)
     outvar[cc]=total(mk_finite(slice_var[pix]*(1./slice_var[pix])^2)*mask)

  endfor

  ;;norm 
  outspc=outspc/outwgt
  outvar=outvar/(outwgt)^2.

  if keyword_set(write) then begin
     
     mwrfits, outspc, write, /cre
     mwrfits, sqrt(outvar), write
     mwrfits, wave, write
     x_specplot, write, inflg=2
  endif
  
  if keyword_set(show) then x_specplot, outspc, wav=wave, sig=sqrt(outvar)

end









