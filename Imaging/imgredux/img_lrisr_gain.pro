;+
;
; Apply the gain in LRISr images.
;
;
;-


pro img_lrisr_gain, adu_image, electron_image, outgain=outgain, outscale=outscale, side=side, notweakgain=notweakgain

  common ccdinfo, xfull,yfull,nchip,namp,biaslev,satur, goodata
  
  ;;define gain value as from lris page
  if(side eq 'R') then begin
      ;;gain value 
      gain_value=[1.98,2.17]
  endif else stop
  

;;tweak the gain to level the two ampli in each chip
;;--------------------------------------------------
  
  med_amp=fltarr(namp)
  
  for amp=0, namp-1 do begin
      xs=1024*amp
      xe=1024*(amp+1)-1

      ;;consider only data regions [430:3750,740:3160] 
      if(xs LT goodata[0]) then xs=goodata[0]
      if(xe GT goodata[1]) then xe=goodata[1]
      ;;try to mask bright stuff
      djs_iterstat, adu_image[xs+100:xe-100,goodata[2]:goodata[3]]*gain_value[amp], median=medaa
      med_amp[amp]=medaa
  endfor
  
  
  if ~keyword_set(notweakgain) then begin
     ;;correct
     splog, 'Tweak gain... Corr_12 :', med_amp[0]/med_amp[1]
     
     ;;check if resonable 
     if(med_amp[0]/med_amp[1] gt 1.2 or med_amp[0]/med_amp[1] lt 0.8) then begin
        splog, 'Variation for chips ', 1, ' bigger than 20%!'
        splog, 'Restore gain ', gain_value[1]
        med_amp[0]=1.
        med_amp[1]=1.
     endif
     
     gain_value[1]=gain_value[1]*med_amp[0]/med_amp[1]
     
  endif
  
  ;;apply final gain
  electron_image=adu_image
  for amp=0, namp-1 do begin
      xs=1024*amp
      xe=1024*(amp+1)-1
      electron_image[xs:xe,0:yfull-1]=adu_image[xs:xe,0:yfull-1]*gain_value[amp]
  endfor
    
  ;;set out gain/scale
  outgain=gain_value
  outscale=[1.,1.]
  
end
