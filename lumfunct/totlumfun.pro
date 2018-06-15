;+
;
; Take as input a file with luminosity function (see Reddy2009_z2)
; and outputs the total number of galaxies brighter than a given 
; absolute magnitude. 
;
; evalmag  --> AB absolute magnitude 
; type     --> any luminosity function with corresponding file type.dat
;
; Output  number of galaxies > evalmag (h70^3 Mpc^-3)
;
; Simple linear interpolation in log log
; Integration is done in linear space
;-





function evallf, x

common lf, absmag, phi

  y=interpol(alog10(phi),absmag,x)
  return, 10^y

end



function totlumfun, evalmag, type

common lf, absmag, phi

;;read 
path=getenv("MIDL")+"/lumfunct/"
readcol, path+type+".dat", minmag, maxmag, phi1e3, error, /silent

;;assign common
absmag=(maxmag-minmag)*0.5+minmag
phi=phi1e3*1D-3


;;integrate
out=qromb('evallf',-28,evalmag) 

return, out

end
