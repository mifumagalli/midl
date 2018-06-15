



pro img_imacs_gain, adu_image, electron_image, header, outgain=outgain, $
                   outscale=outscale, side=side


  ;;grab nominal gain
  egain=sxpar(header,"EGAIN")
  
  
  ;;tweak the gain to level in each
  ;;chip. Consider central chip (number 2) the reference.
  ;;Note that this is a somewhat arbitrary initialization as the gain
  ;;of mosaic 3 are not published. WTF. Alos, note that the order of
  ;;this array is relative to the geometry of the ccd as in img_oscan
  
  
  ;;gain_value=[0.924153,0.958004,0.870000,0.860774,0.89000,0.860380,0.896955,0.933331]    ;;V
  ;;gain_value=[0.984153,0.958004,0.870000,0.900774,0.942276,0.866380,0.896955,0.953331]   ;;R&B
  gain_value=[0.944153,0.948004,0.870000,0.900774,0.902276,0.856380,0.896955,0.953331]   ;;I

  med_amp=fltarr(8)
  outgain=fltarr(8)
  
  for amp=0, 7 do begin
     ;;set the y
     ymin=0
     ymax=4095
     ;;find position 
     if(amp le 3) then begin
        ys=ymin
        ye=ymax
        xs=8192-(amp+1)*2048
        xe=8191-amp*2048
     endif else begin 
        ys=ymin+ymax+1
        ye=ymax*2+1
        xs=8192-(amp-4+1)*2048
        xe=8191-(amp-4)*2048
     endelse
     med_amp[amp]=1.;djs_median(adu_image[xs:xe,ys:ye])*gain_value[amp]
  endfor
  
  for amp=0, 7 do begin
     
     ;;correct (reference central chip)
     outgain[amp]=gain_value[amp]*med_amp[2]/med_amp[amp]
     splog, 'Tweak gain... ', amp, outgain[amp], ' Corr: ', med_amp[amp]/med_amp[2]
     
  endfor

  ;;check if resonable 
  nogood=where(outgain/gain_value gt 1.2 or  outgain/gain_value lt 0.8, ncheck)
  if(ncheck gt 0) then begin
     splog, 'Variation for chips ', nogood+1, ' bigger than 20%!.'
     splog, 'Restore gain ', gain_value[nogood]
     outgain[nogood]=gain_value[nogood]
  endif
  
  ;;apply gain
  electron_image=adu_image
  
  for amp=0, 7 do begin
     ;;set the y
     ymin=0
     ymax=4095
     ;;find position 
     if(amp le 3) then begin
        ys=ymin
        ye=ymax
        xs=8192-(amp+1)*2048
        xe=8191-amp*2048
     endif else begin 
        ys=ymin+ymax+1
        ye=ymax*2+1
        xs=8192-(amp-4+1)*2048
        xe=8191-(amp-4)*2048
     endelse
     electron_image[xs:xe,ys:ye]=adu_image[xs:xe,ys:ye]*outgain[amp]

  endfor
  
  ;;for now lbc has scale 1
  outscale=outgain-outgain+1

end