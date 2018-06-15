;+
;
; Read a spectrum written by class for the IRAM 30 mt
;
; freq in GHz
; temp in K
; vel  in km/s
; lambda in mm
;
; plot use spec_plot
;
;-




pro readiram, file, freq=freq, temp=temp, vel=vel, $
              lambda=lambda, plot=plot

  
  temp=mrdfits(file,0,head,/fsc,/sil)
  
  ;;calibrate in frequency
  RESTFREQ=fxpar(head,"RESTFREQ")
  CRVAL1=fxpar(head,"CRVAL1")
  CRPIX1=fxpar(head,"CRPIX1")
  CDELT1=fxpar(head,"CDELT1")
  DELTAV=fxpar(head,"DELTAV")
  VLSR=fxpar(head,"VLSR")
  
  INDX=findgen(n_elements(temp))
  freq=(RESTFREQ+CRVAL1+(INDX-CRPIX1)*CDELT1)/1d9 ;GHz
  vel=VLSR+(INDX-CRPIX1)*DELTAV*1D-3              ;km/s
  lambda=2.99792458000D11/(freq*1d9)              ;mm

  if keyword_set(plot) then begin
     tmp="tmpSDBJfedwvr6iuNLHBJCD.fits"

     mwrfits, temp, tmp, /cre, /sil
     mwrfits, temp-temp, tmp
     mwrfits, lambda, tmp
     
     x_specplot, tmp, inflg=2
     
     spawn, "rm -f "+tmp
     
  endif
  
end
