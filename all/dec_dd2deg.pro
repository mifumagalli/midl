
;convert a number ddmmss.dec into deg DEC


FUNCTION dec_dd2deg, dec

dec=strtrim(dec,2)
ra=replicate('0:00:00.0',n_elements(dec))


x_radec, ra, dec, ra_deg, dec_deg

return, dec_deg

END

  
