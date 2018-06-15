function wilsonscore, p, n, confint=confint, norm=norm, verbose=verbose 

; JT 111510 
; JXP 061512 -- Made some changes to handle bigger numbers.
; calculates binominal confidence intervals (estimates mean of binomial distribution) 
; from observations. 
; p = "hit rate" measured in data, ie. 20 trials, 19 hits = 0.95 
; n = number of trials 
; confint = confidence interval total width. Calculated score will contain 
;  		that fraction of the total probabilty. 95% by default. 



  ; initially, assume 95% interval 
  if (not keyword_set(confint)) then confint = 0.95 ;; 2-sigma equivalent
 
  z = gauss_cvf( (1. - confint) / 2.d ) 
  if (keyword_set(verbose)) then  print, 'Assuming confidence interval', confint, z 

 if (keyword_set(norm)) then begin 
  ; now try the normal approximation 
  hi = p + z * (p * (1.d - p) / n ) ^ 0.5 
  lo = p - z * (p * (1.d - p) / n ) ^ 0.5 
  if (keyword_set(verbose)) then print, 'Calculating Normal Approximation',  lo, p, hi 
 endif else begin 
  ; now do it with Wilson score interval  
  center = (p + 1. / 2.d / n * z^2) / (1. + z^2 / n) 
  hi = (p + 1. / 2.d / n * z^2 + z * ( (p * (1.-p) / n) + (z^2 / 4. / n^2.) )^0.5 ) / (1. + z^2 / n) 
  lo = (p + 1. / 2.d / n * z^2 - z * ( (p * (1.-p) / n) + (z^2 / 4. / n^2.) )^0.5 ) / (1. + z^2 / n) 
  if (keyword_set(verbose)) then print, 'Calculating Score Wilson Interval',  lo, p, hi 
 endelse 

  return, [lo, center, hi] 

  

end 
