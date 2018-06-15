;+
;PURPOSE	
;	to compute a dust correction using milky way dust curve 
;SYNTAX
;	corr=rdust(lambda,a_lam,A_v)
;INPUTS
;	lambda: the wavelength [Ang]
;	A_v: V-band extinction
;
;OUTPUTS
;       a_lam: the dust extinction in magnitude 
;	
;NOTES
;	requires use of XIDL to get the dust curve and the $XIDL_DIR 
;	environment variable to be set
;
;Written by R. da Silva, UCSC, 8-26-09
;Modified by MF, June 2010 
;-
FUNCTION dust_law, lambda, A_v

  if N_PARAMS() NE 2 then begin
     print, 'Synatx - corr=rdust(lambda, flux, A_v)'
     return, -1
  endif  
    
  ;; Calculating Correction
  dir=getenv('XIDL_DIR')+'/Dust/' 
  file=dir+'MW_dust.dat'        ;grabbing dust curve
  readcol, file, lam, xi,/silent ;and reading it in
  
  xis=interpol(xi, lam, lambda)	;interpolated to lambda for evaluation
  
  a_lam=xis*A_v			;converting xi to A_lambda
  
  ;;correction=10.^(.4*a_lam)	;finding the correction to apply to flux
  ;;                              
  ;;recall A_lambda=-2.5*alog10(f_observed/f_intrinsic)
  ;;corr=flux*correction		;applying the correction
  ;;  return, corr			;returning the values

  return, a_lam

end

