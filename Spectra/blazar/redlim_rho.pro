;+
;
;
;  Compute the log of rho_0 following equation 4 in Sbarufatti+2006
;
;  redshift   the redshift array 
;  mrnuc      the nuclear magnitude
;  mstarz     the host galaxy magnitude array
;
;  
;-


pro redlim_rho, redshift, mrnuc, mstarz, logrho=logrho

  ;;distance
  distance_calculator, redshift, lumdist, /lum, /wmap7, /sil 
  
  ;;k-correct as in Wisotzki+00
  kcorr=redshift-redshift
  
  low=where(redshift le 0.6)
  kcorr[low]=-0.07-2.5*(1+0.37)*alog10(1+redshift[low])
  
  if(max(redshift) gt 0.6) then begin
     high=where(redshift gt 0.6)
     kcorr[high]=-0.48-2.5*(1-0.60)*alog10(1+redshift[high])
  endif

  mnucz=(mrnuc+5-5*(alog10(lumdist)+6)+kcorr)
  logrho=-0.4*(mnucz-mstarz)

end
