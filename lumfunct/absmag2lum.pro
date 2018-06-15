;+
;
;Procedure that converts absolute AB mag into a luminosity or other
;way around.
;
;
;LUM ---> if set, input is assumed to be a MAG (AB) and output is LUM
;         (erg/s/Hz)
;
;MAG ---> if set, input is assumed to be a LUM (erg/s/Hz) and output is
;         MAG (AB)
;
;-

pro absmag2lum, input, output, LUM=lum, MAG=mag


;;make everything a double
input=1D0*input


;Absolute magnitude is observed magnitude at D_10=10 pc.
;M=-2.5 Log (F/4 pi D_10^2)-48.6
;


;;get absolute magnitude from luminosity
if keyword_set(mag) then output=-2.5*alog10(input)+51.595006

;;get luminosity from absolute magnitude
if keyword_set(lum) then output=10^(-0.4*(input-51.595006))

end
