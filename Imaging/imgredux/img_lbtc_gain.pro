



pro img_lbtc_gain, adu_image, electron_image, outgain=outgain, $
                   outscale=outscale, side=side


;;tweak the gain to level in each
;;chip. Consider central chip (number 2,index1) the reference.
if(side eq 'R') then gain_value=[2.13,2.14,2.08,2.09]
if(side eq 'B') then gain_value=[2.06,2.09,1.96,1.98]


med_amp=fltarr(4)
  
;;get the scale for each ccd
med_amp[0]=djs_median(adu_image[0:2047,0:4607])*gain_value[0]
med_amp[1]=djs_median(adu_image[2122:4169,0:4607])*gain_value[1]
med_amp[2]=djs_median(adu_image[4244:6291,0:4607])*gain_value[2]
med_amp[3]=djs_median(adu_image[770:5377,4683:6730])*gain_value[3]


;;correct (reference central chip)
splog, 'Tweak gain... Corr_12: ', med_amp[0]/med_amp[1]
splog, 'Tweak gain... Corr_32: ', med_amp[2]/med_amp[1]
splog, 'Tweak gain... Corr_42: ', med_amp[3]/med_amp[1]

outgain=fltarr(4)
outgain[0]=gain_value[0]*med_amp[1]/med_amp[0]
outgain[1]=gain_value[1]
outgain[2]=gain_value[2]*med_amp[1]/med_amp[2]
outgain[3]=gain_value[3]*med_amp[1]/med_amp[3]

;;check if resonable 
nogood=where(outgain/gain_value gt 1.2 or  outgain/gain_value lt 0.8, ncheck)
if(ncheck gt 0) then begin
   splog, 'Variation for chips ', nogood+1, ' bigger than 20%!.'
   splog, 'Restore gain ', gain_value[nogood]
   outgain[nogood]=gain_value[nogood]
endif

;;apply gain
electron_image=adu_image
electron_image[0:2047,0:4607]=adu_image[0:2047,0:4607]*outgain[0]
electron_image[2122:4169,0:4607]=adu_image[2122:4169,0:4607]*outgain[1]
electron_image[4244:6291,0:4607]=adu_image[4244:6291,0:4607]*outgain[2]
electron_image[770:5377,4683:6730]=adu_image[770:5377,4683:6730]*outgain[3]

;;for now lbc has scale 1
outscale=outgain-outgain+1

end
