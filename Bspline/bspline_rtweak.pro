function bspline_rtweak, image, ra, dec, rdc_vec, $
 invvar=invvar, rval=rval, ntweak=ntweak, _Extra=kw_for_radial

; WRITTEN: abolton@cfa

if (not keyword_set(rval)) then rval = .5
if (not keyword_set(ntweak)) then ntweak = 5
if (n_elements(invvar) ne n_elements(image)) then invvar = 1. + 0.*image

tval = (!dpi * dindgen(360) / 180.d0) - !dpi
rr = rval # replicate(1., n_elements(tval))
tt = replicate(1., n_elements(rval)) # tval
xx = rr * cos(tt)
yy = rr * sin(tt)

ra_cen = rdc_vec[0]
dec_cen = rdc_vec[1]
for i = 0, ntweak-1 do begin
  delta_ra = ra - ra_cen
  delta_dec = dec - dec_cen
  r = sqrt(delta_ra^2 + delta_dec^2)
  theta = atan(delta_dec, delta_ra)
  rset = bspline_radial(r, theta, image, invvar=invvar, yfit=yfit, $
   _Extra=kw_for_radial)
  ff = bspline_radial_valu(rr, tt, rset)
  xff = total(xx*ff, 2) / total(ff, 2)
  yff = total(yy*ff, 2) / total(ff, 2)
  ra_cen = ra_cen + mean(xff)
  dec_cen = dec_cen + mean(yff)
endfor

return, [ra_cen, dec_cen]
end

