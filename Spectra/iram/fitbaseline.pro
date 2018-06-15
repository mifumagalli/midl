;+
;
;   Fit a polinomial baseline to a spectrum, skipping a window
;
;-

pro fitbaseline, velo, spec, deg=deg, base=base, $
                   window=window
  
  ;;suppress fit in window
  error=velo-velo+1
  error[where(velo gt window[0] and velo lt window[1])]=1d30

  result = POLY_FIT(velo,spec,deg,yfit=base,MEASURE_ERRORS=error)
  
end
