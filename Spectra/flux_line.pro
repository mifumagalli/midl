;+
;
;  Take a 1D spectrum, and min max range
;  and do a simple estimate of the line flux and EW.
;  This can be improved considering the error.
;
;  inflg      the inflag used to reading the spectrum
;  lamobs     observed lambda of transition
;  lamres     rest frame lambda of transition
;  lamrange   the lambda interval used to evaluate the line 
;  sky        an array of lambda gap, and width used to estimate 
;             consitnuum 
;  
;
;-




pro flux_line, spec, inflg=inflg, lamobs=lamobs, lamres=lamres, $
               lamrange=lamrange, sky=sky


  if not keyword_set(inflg)  then inflg=2
  if not keyword_set(lamobs) then lamobs=5000.
  if not keyword_set(lamres) then lamres=5000.
  if not keyword_set(lamrange) then lamrange=[lamobs-10,lamobs+10]
  if not keyword_set(sky) then sky=[10.,50.]

  
    
  ;; redshift and D luminosity
  red=lamobs/lamres-1
  distance_calculator, red, dlum, /lum, /wmap5, /sil
  
  splog, 'Redshift ', red
  splog, 'Distance (Mpc) ',dlum
  
  
  ;;read spec
  flux=x_readspec(spec,inflg=inflg,head=head,sig=sig,wav=wav)

  
  ;;find crude sky level
  skyright_min=lamrange[1]+sky[0]
  skyright_max=skyright_min+sky[1]
  skyleft_max=lamrange[0]-sky[0]
  skyleft_min=skyleft_max-sky[1]

  
  skylev_left=median(flux[where(wav gt skyleft_min and wav lt skyleft_max)]) 
  skylev_right=median(flux[where(wav gt skyright_min and wav lt skyright_max)]) 

  skylev=0.5*(skylev_left+skylev_right)

  splog, 'Sky left-right ',  skylev_left, skylev_right
  splog, 'Sky level ', skylev


  ;;find line properties
  inline=where(wav gt lamrange[0] and wav lt lamrange[1])
  dlambda=wav-shift(wav,1)
  ;;integrate over line and subtract sky
  tot_flux=total((flux[inline]>0-skylev)*dlambda[inline])
  tot_lum=tot_flux*(4*!pi*(dlum*1D6*3.0856802D18)^2)
  
  splog, 'Total line flux ', tot_flux
  splog, 'Total line lum  ', tot_lum


  ;;find EW
  ew_obs=total((1-(flux[inline]>0)/skylev)*dlambda[inline])
  ew_res=ew_obs/(1+red)
  splog, 'Obs EW ', ew_obs 
  splog, 'Res EW ', ew_res
  

  




  ;;plot
  maxx=lamrange[1]+sky[0]+sky[1]+5
  minx=lamrange[0]-sky[0]-sky[1]-5

  ;;spec
  plot, wav, flux, xrange=[minx,maxx], psym=10, $
          xtitle='lambda', ytitle='flux'
  oplot, wav, sig, psym=10, color=250

  ;;line
  oplot, [lamobs,lamobs], [-1D50,1D50], line=2
  oplot, lamrange, [skylev,skylev], line=1, color=50, thick=8
  
  ;;sky
  oplot, [skyleft_min,skyleft_max], [skylev,skylev], line=1, color=150, thick=8
  oplot, [skyright_min,skyright_max], [skylev,skylev], line=1, color=150, thick=8

  


end
