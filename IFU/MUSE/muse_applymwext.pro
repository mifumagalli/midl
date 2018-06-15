;+
;
; Apply dereddeing correction for the MW absorption
;
; datacube -> reconstructed MUSE datacube
; mwext -> File with the amount of extension as a function of wavelenght.
; outfile -> name of the output file
;
;-

pro muse_applymwext, datacube, mwext, outfile

prihdr = headfits(datacube, exten=0, /silent)
data = mrdfits(datacube, 1, head1, /silent)
noise = mrdfits(datacube, 2, head2, /silent)

wave = fxpar(head1, 'CRVAL3')+indgen(fxpar(head1, 'NAXIS3'))*fxpar(head1, 'CD3_3')
readcol, mwext, lambda_dust, kappa, corr_dust

dustcorr = interpol(corr_dust, lambda_dust, wave , /spline)

Ncol = (size(data))[1]
Nrow = (size(data))[2]

for col=0, Ncol-1 do begin
 for row=0, Nrow-1 do begin
  data[col,row,*] *= dustcorr
 endfor
endfor

mwrfits, 0, outfile, prihdr, /create
mwrfits, data, outfile, head1
mwrfits, noise, outfile, head2

end
