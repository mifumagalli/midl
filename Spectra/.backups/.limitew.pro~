;+
;
;  Given the variance array, a velocity window and the  
;  rest vawe, find the 2sigma upper limits to the line 
;  
;
;  errspec    the spectrum error for normalized data
;  zabs       redshift of absorption
;  velocity   [minv,maxv] to be used for the integration
;  lambda     rest frame velocity of the transition 
;  data       the actual normalized flux
;-

pro limitew, errspec, zabs=zabs, velocity=velocity, $
             lambda=lambda, inflg=inflg, data=data
  
  ;;open the err spec
  err=x_readspec(errspec,0,INFLG=inflg,wav=wav)
  if keyword_set(data) then flux=x_readspec(data,0,INFLG=inflg,wav=wavflux)
  
  
  ;;find the relevant region for the fit
  lobs=(1+zabs)*lambda
  wvel=(wav-lobs)*299792.458/lobs
  vrange=where(wvel le velocity[1] and wvel ge velocity[0])
  
  ;;measure the rest frame EW
  linstr=x_setline(lambda)
  ew=x_calcew(wav,1-2*err,minmax(vrange),/FPIX)/(1+zabs)
  splog, "Rest frame 2 sigma EW ", ew

  ;;get column density with direct integration of opacity 
  minl=lobs*velocity[0]/299792.458+lobs
  maxl=lobs*velocity[1]/299792.458+lobs
  lrange=where(wav ge minl and wav le maxl)

  x_aodm, wav[lrange], 1-2.*err[lrange], lrange-lrange, lambda, colm, sig_colm, /LOG  
  splog, "Column density limit ", colm

  ;;plot 
  prange=where(wvel le 3*velocity[1] and wvel ge 3*velocity[0])
  
  plot, wvel[prange], 1-2*err[prange], psym=10, yrange=[-0.2,1.2]
  oplot, [-1d5,1d5], [1.,1.], color=fsc_color("green"), line=1 
  oplot, [-1d5,1d5], [0.,0.], color=fsc_color("green"), line=1
  oplot, [0,0], [-100.,100.], color=fsc_color("red"), line=1
  
  oplot, (wav[lrange]-lobs)*299792.458/lobs, 1-2*err[lrange], psym=5, color=250

  if keyword_set(data) then oplot, wvel[prange], flux[prange], $
                                   color=fsc_color("blue"), psym=10
  

end
