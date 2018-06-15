;+
;
; Build the frequency axis for a spectrum, given the header
;
; freq in GHz
; 
; multi     set this is passing a structure instead of single spectrum
;-




pro iramfreq, head, freq=freq, multi=multi
  
  ;;calibrate in frequency
  RESTFREQ=fxpar(head,"RESTFREQ")
  CRVAL1=fxpar(head,"CRVAL1")
  CRPIX1=fxpar(head,"CRPIX1")
  CDELT1=fxpar(head,"CDELT1")
  if keyword_set(multi) then NCHANN=fxpar(head,"MAXIS1")  else $
  NCHANN=fxpar(head,"NAXIS1")
  
  INDX=findgen(NCHANN)
  freq=(RESTFREQ+CRVAL1+(INDX-CRPIX1)*CDELT1)/1d9 ;GHz
   
end
