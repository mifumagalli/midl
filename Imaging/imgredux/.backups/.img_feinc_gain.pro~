;+
;
; Apply the gain in ESI images.
;
; Gain (low) 	
; 1.29 e-/DN
; Gain (high) 	
; 0.5 e-/DN
;
;-


pro img_esi_gain, adu_image, electron_image, header, outgain=outgain, outscale=outscale

  common ccdinfo, xfull,yfull,nchip,namp,biaslev,satur, goodata
  
  
  ;;set gain level
  speed=sxpar(header,"CCDGAIN")
  if(rstring(speed) eq "low") then gain_value=1.29 else gain_value=0.5
  
  ;;apply final gain
  electron_image=adu_image*gain_value
  
  ;;set out gain/scale
  outgain=gain_value
  outscale=1.
  
end
