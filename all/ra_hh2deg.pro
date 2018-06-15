;convert a number hhmmss.dec into deg RA


FUNCTION ra_hh2deg, ra

ra=strtrim(ra,2)
dec=replicate('0:00:00.0',n_elements(ra))
x_radec, ra, dec, ra_deg, dec_deg

return, ra_deg

END

  
