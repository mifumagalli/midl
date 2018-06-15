;+
;
;
;  This is a tiny script that makes the atmospheric 
;  extinction as a function of lambda from the AM terms
;
;
;-




pro lbt_ext


;;read AM extinction
readcol, 'LBTC_Extinction.dat', name, mag_ext, mag_ext_err, lambda_nm, $
  format='A,F,F,F'

;;make AA
lambda_aa=lambda_nm*10.
ext=10^(0.4*mag_ext)


forprint, lambda_aa, ext, textout='LBTC_atmextinction.dat'



end
