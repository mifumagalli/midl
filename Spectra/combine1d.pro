;+
;
;
; Combine a set of spectra using the x_combspec procedure
; Make additional tests to make sure the wavelenghts are fine
;
;
;  spect      list of spectra to coadd
;  outname    output coadded spectrum
;  nrmflux    set to window range to use in flux normalization
;-


pro combine1d, spect, outname, nrmflux=nrmflux, inflg=inflg
  
  if ~keyword_set(nrmflux) then nrmflux=[5500.,7500.]


  ;;load first 
  if keyword_set(inflg) then begin
     flux_one=x_readspec(spect[0]+'.fits',inflg=inflg,$
                         sig=error_one, wav=wav_one, head=head) 
  endif else begin
     flux_one=x_readspec(spect[0]+'_F.fits',fil_sig=spect[0]+'_E.fits',$
                         sig=error_one, wav=wav_one, head=head) 
  endelse
  
  
  ;;prepare to stack
  flux_all=fltarr(n_elements(flux_one),n_elements(spect))
  var_all=fltarr(n_elements(flux_one),n_elements(spect))

  
  ;;add first
  flux_all[*,0]=flux_one
  var_all[*,0]=error_one^2
  
  ;;load and stack next 
  for ss=1, n_elements(spect)-1 do begin
     
     if keyword_set(inflg) then begin
        flux_two=x_readspec(spect[ss]+'.fits',inflg=inflg,$
                            sig=error_two, wav=wav_two, head=head2) 
     endif else begin
        flux_two=x_readspec(spect[ss]+'_F.fits',fil_sig=spect[ss]+'_E.fits',$
                            sig=error_two, wav=wav_two, head=head2) 
     endelse
     

     date=fxpar(head2,'UT-DATE')
     
     splog, "Adding ", spect[ss], ' ', date
     splog, "Wavelenght error ", total(wav_one-wav_two)
     splog, "Elements in array ",  n_elements(wav_one), n_elements(wav_two)
    
     ;;if there is error in wavelength, stop
     if(total(wav_one-wav_two) ne 0) then begin
     
        splog, 'Check interpolation '
        
        plot, wav_two, flux_two, xrange=[3500,3600]
   
        ;;do linear interpolation on first array
        flux_two=interpol(flux_two,wav_two,wav_one)
        error_two=interpol(error_two,wav_two,wav_one)
        wav_two=wav_one
        
        oplot, wav_two, flux_two, color=fsc_color('red'), line=1
        
        stop
        

     endif
     
     ;;if there is no error in wave, but it's just a truncation issue
     ;;cut the tail
     if(n_elements(wav_one) ne n_elements(wav_two)) then begin
       
        ;;wav_one is smaller: cut tail
        if(n_elements(wav_one) lt n_elements(wav_two)) then begin
           splog, 'Dump the tail'
           flux_two=flux_two[0:n_elements(wav_one)-1]
           error_two=error_two[0:n_elements(wav_one)-1]
        endif else begin
           splog, 'Fill remaining with 0s'
           tmpf_two=flux_two
           tmpe_two=error_two
           flux_two=flux_one-flux_one
           error_two=error_one-error_one
           flux_two[0:n_elements(wav_two)-1]=tmpf_two
           error_two[0:n_elements(wav_two)-1]=tmpe_two 
        endelse
        
     endif
     
     ;;add to stack
     flux_all[*,ss]=flux_two
     var_all[*,ss]=error_two^2
     
     ;;update header
     fxaddpar, head, 'COMMENT', 'Added 1D spectrum from '+date  
	
  endfor

  ;;prepare final stack
  x_combspec, flux_all, var_all, fflux, fvar, wave=wav_one, nrmflux=nrmflux
  

  splog, 'Check flux_one/flux_all '
  plot, wav_one, fflux/flux_one, xrange=[3100,8000]
  stop
  
  splog, 'Check error_one/error_all '
  plot, wav_one, error_one/sqrt(fvar), xrange=[3100,8000] 
  stop

  splog, 'Check SN gain'
  plot, wav_one, fflux/sqrt(fvar), xrange=[3100,8000] 
  oplot, wav_one, flux_one/error_one, color=fsc_color('red')
  stop

  ;plot, wav_one, flux_one, xrange=[3100,8000]
  ;oplot, wav_two, flux_two, color=250


  ;;write output
  if ~keyword_set(inflg) then begin
     mwrfits, fflux, outname+'_F.fits', head, /create
     mwrfits, sqrt(fvar), outname+'_E.fits', head, /create
     ;;check
     x_specplot, outname+'_F.fits', /auto, /bl
     x_specplot, spect[0]+'_F.fits', /auto, /bl
  endif else begin

     if(inflg eq 2) then begin
        mwrfits, fflux, outname+'.fits', head, /create
        mwrfits, sqrt(fvar), outname+'.fits'
        mwrfits, wav_one, outname+'.fits'
        
        ;;check
        x_specplot, outname+'.fits', /bl, inflg=inflg
  
     endif else stop


  endelse

 

end
