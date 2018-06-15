;+
;
; Take a low resolution spectrum (inflg=2) and the constinuum
; (derive e.g. using x_continuum) and produce the normalize spectrum
;
;
;-




pro normlowres, spec, cont, out


  ;;read in
  flux=mrdfits(spec,0,hd0)
  error=mrdfits(spec,1,hd1)
  wave=mrdfits(spec,2,hd2)
  contlev=mrdfits(cont,0)

  ;;normalize
  norm_flux=flux/contlev
  norm_error=error/contlev


  ;;write 
  mwrfits, norm_flux, out, hd0, /create
  mwrfits, norm_error, out, hd1
  mwrfits, wave, out, hd2
  



end
