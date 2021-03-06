;+
;
;
;  Do the actual combine, after registering to the same frequency 
;  frame, applying heliocentric correction  
;  
;
;-


pro scancoadd_add, filein, out=out

  ;;loop and coadd
  for i=0, n_elements(filein)-1 do begin

     str=mrdfits(filein[i]+'_proc.fits',1)
     nscan=n_elements(str)
     nchan=n_elements(str[0].spectrum)

     
     ;;loop over scans

     for ss=0, nscan-1 do begin
        
        ;;initialize running stack
        if (i eq 0 and ss eq 0) then begin
           stack_velocity=str[0].VELCORR
           stack_weight=replicate(str[0].OBSTIME/sqrt(str[0].TSYS),nchan)
           stack_spectrum=str[0].SPECSUB*stack_weight
           stack_wave=str[0].wavecor

           ;;find channel separation
           deltav=abs(stack_velocity-shift(stack_velocity,1))
           deltav=median(deltav[1:nchan-2])
           print, 'Velocity channel ', deltav
        endif
        

        ;;do the coadd with an ineger shift
        relativev=median(stack_velocity-str[ss].VELCORR)
        mumrel=round(relativev/deltav)
        shiftvel=median(stack_velocity-shift(str[ss].VELCORR,mumrel))
        if(abs(shiftvel) gt 0.5*deltav) then begin
           splog, 'Shift exceed tolerance!'
           stop
        endif
        
        ;;update stack
        weight=replicate(str[ss].OBSTIME/sqrt(str[ss].TSYS),nchan)
        ;;mask the region for shift
        if(mumrel gt 0) then weight[0:mumrel]=0
        if(mumrel lt 0) then weight[nchan+mumrel:nchan-1]=0
        stack_spectrum+=shift(str[ss].SPECSUB,mumrel)*weight
        stack_weight+=weight
     endfor
     
  endfor

  
  ;;do the final stack
  stack_spectrum=stack_spectrum/stack_weight
  

  ;;write
  mwrfits, stack_spectrum, out, /crea
  mwrfits, replicate(sqrt(mean(stack_spectrum^2)),nchan), out
  mwrfits, stack_wave, out
  mwrfits, stack_velocity, out
  
end
