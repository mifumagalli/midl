;+
;
;  Load the muse datacube and error, parsing headers 
;
;  filename - name of the cube 
;  wave     - wave in A (heliocentric corrected and in vacuum)
;  flux     - flux cube (1d-20 erg/s/cm^2/A)
;  variance - error cube (1d-20 erg/s/cm^2/A)^2 
;  unwrap   - run the unwrap procedure 
;  astro    - the astrometry structure 
;  mmlam    - min max lambda for wrapping 
;  obswave  - original wave (in air, no helio corr)
;-

pro muse_loadcube, filename, wave=wave, flux=flux, var=var, path=path,$
                   unwrap=unwrap, astro=astro, mmlam=mmlam, obswave=obswave


  if ~keyword_set(path) then path=""

  ;;load data 
  null=mrdfits(path+filename,0,mainheader)
  flux=mrdfits(path+filename,1,sciheader)
  var=mrdfits(path+filename,2,varheader)

  ;;load astro 
  extast, sciheader, astro

  ;;reconstruct wavelength 
  sz=size(flux)
  delta_lambda=fxpar(sciheader,"CD3_3")
  zero_lambda=fxpar(sciheader,"CRVAL3")
   
  ;;this is in air 
  wave=make_array(sz[3],/index)*delta_lambda+zero_lambda
  
  ;;check for helio correction 
  heliocorr=fxpar(sciheader,"HELIO")
  
  if (heliocorr ne 0) then begin
     
     ;;helio correction is applied.. undo
     obswave=wave/heliocorr
     
     ;;;;go to vacuum
     ;;;;airtovac, wave
     
  endif else begin
     
     splog, 'HELIO CORR NOT APPLIED!!!!'
     obswave=wave ;;this is air no helio
     
  endelse 
  
  ;;unwrap 
  if keyword_set(unwrap) then begin

     splog, 'Unwrapping the cube'
     
     uwrapname=strmid(filename,0,strpos(filename,".fits"))

     ifu_unwrap, wave, flux, twodflux, lambda=mmlam
     mwrfits, twodflux, uwrapname+"_sciunwr.fits", /crea
     undefine, twodflux
     
     ifu_unwrap, wave, var, twodvar, lambda=mmlam
     mwrfits, twodvar, uwrapname+"_varunwr.fits", /crea
     undefine, twodvar

  endif
  


end
