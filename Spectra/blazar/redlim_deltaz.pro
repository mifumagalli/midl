;+
;
;  Compute the delta_z function from equation 3 of 
;  Sbarufatti et al. 2006
;
;  alpha         the blazar power law
;  delta_z       output
;  lambda_zero   wavelength of normalization
;  lambda_line   wavelength of transiton
;
;-



pro redlim_deltaz, wav, flux, redshift, alpha=alpha, delta_z=delta_z, $
                   lambda_zero=lambda_zero, lambda_line=lambda_line
  
  if ~keyword_set(alpha) then begin

     fit=linfit(alog10(wav[where(wav gt 3500 and wav lt 5500)]), $
                alog10(flux[where(wav gt 3500 and wav lt 5500)]),YFIT=yfit)
     
     plot, wav, flux, xrange=[3500,5500.], xstyle=1
     oplot, wav, 10^(alog10(wav)*fit[1]+fit[0]), line=1, $
            color=fsc_color("red"), thick=5.
     
     alpha=fit[1]
     splog, "Spectral index ", alpha
     
  endif

  ;;set nucleus normalization
  nucleus=wav^alpha
  delta_z=redshift-redshift

  ;;load early type template 
  readcol, getenv('MIDL')+'/Spectra/blazar/elliptical_template.ascii',$
           gal_wav, gal_flux
  
  ;;loop over redshift 
  for zz=0, n_elements(redshift)-1 do begin
     
     ;;find the rho at each redshift
     gal_wav_z=gal_wav*(1+redshift[zz])
     gal_flux_int=interpol(gal_flux,gal_wav_z,wav)

     ;;get rho_lambda
     rho_lambda=(nucleus/gal_flux_int)
     norm=interpol(rho_lambda,wav,lambda_zero)
     rho_lambda=rho_lambda/norm[0]
     
     ;;evaluate at CaII
     this_lambda=lambda_line*(1+redshift[zz])
     delta_z[zz]=interpol(rho_lambda,wav,this_lambda)
     
  endfor
  

end
