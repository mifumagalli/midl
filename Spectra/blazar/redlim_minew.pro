;+
;
;   Find the minimum EW in a spectrum following the method by 
;   Sbarufatti 2005
;
;   wav         wavelength array 
;   spectrum    flux
;   minew       in output the min EW
;   regions     a input file with min/max wave to examine (first line) and 
;               spectral regions to exclude
;
;-


pro  redlim_minew, wav, spectrum, minew=minew, regions=regions


  ;;use bin of 20 A
  delta=20. 

  ;;readinfo
  readcol, regions, min_l, max_l
  nbin=(max_l[0]-min_l[0])/delta

  splog, "Compute minimum EW of nbin ", nbin
  minew_arr=fltarr(nbin-2)
  
  ;;contruct pixel mask
  mask=wav-wav
  dwav=wav-shift(wav,1.)
  
  ;;loop over mask
  for mm=1, n_elements(max_l)-1 do begin

     ;;find pixels
     center=where(wav gt min_l[mm] and wav le max_l[mm])
     mask[center]=1.
 
  endfor

  ;;loop over bins
  for ee=0, nbin-3 do begin

     ;;find pixels
     left=where(wav gt min_l[0]+ee*delta and wav le min_l[0]+(ee+1)*delta)
     right=where(wav gt min_l[0]+(ee+2)*delta and wav le min_l[0]+(ee+3)*delta)
     center=where(wav gt min_l[0]+(ee+1)*delta and wav le min_l[0]+(ee+2)*delta)
     
     ;;skip if bad
     if total(mask[center]) gt 0. then minew_arr[ee]=0. else begin
        
        ;;find continuum 
        flux_left=spectrum[left]
        flux_right=spectrum[right]
        flux_center=spectrum[center]
        
        lam_left=min_l[0]+(ee+0.5)*delta
        lam_center=min_l[0]+(ee+1.5)*delta
        lam_right=min_l[0]+(ee+2.5)*delta
        
        ;;check continuum 
        px_left_cnt=where(mask[left] lt 1.,npix_left)
        px_right_cnt=where(mask[right] lt 1.,npix_right)

        if(npix_left eq 0 or npix_right eq 0) then minew_arr[ee]=0. else begin
           
           cont_left=median(flux_left[px_left_cnt])
           cont_right=median(flux_right[px_right_cnt])
           cont_center=((cont_right-cont_left)/(lam_right-lam_left))*$
                       (lam_center-lam_left)+cont_left
           
           ;;compute EW 
           minew_arr[ee]=total((flux_center-cont_center)/cont_center*dwav[center])
           
          ; ;;plot
          ; plot, wav, spectrum, xrange=[min_l[0]+ee*delta-2.,$
          ;                              min_l[0]+(ee+3)*delta+2.]
          ; oplot, [lam_left,lam_center,lam_right], $
          ;        [cont_left,cont_center,cont_right], psym=1, symsiz=3., $
          ;        color=fsc_color("red")
           
        endelse
        
     endelse

  endfor

  ;;exclude bad parts
  ew_nozero=minew_arr[where(minew_arr ne 0.)]
  minew=2.*sqrt(total(ew_nozero^2.)/n_elements(ew_nozero))
  
  splog, "Min EW in A ", minew

end
