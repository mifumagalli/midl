;+
;
;function that returns the Pog mag for a give ashin mag
;as described on sdss website
;
;Pogson
;    mag = -2.5 * log10(f/f0)
;    error(mag) = 2.5 / ln(10) * error(counts) / counts
;    To get the error on the counts, see the note on computing count errors below.
;   
;asinh
;    mag = -(2.5/ln(10))*[asinh((f/f0)/2b)+ln(b)]
;    error(mag) = 2.5 / ln(10) * error(counts)/exptime * 1/2b *
;            100.4*(aa + kk * airmass) / sqrt(1 + [(f/f0)/2b]2),
;    where b is the softening parameter for the photometric band 
;    in question and is given in the table of b coefficients below 
;    (for details on the asinh magnitudes, see the paper by Lupton, $
;     Gunn, and Szalay 1999 [AJ 118, 1406]).
;
;
;          asinh Softening Parameters (b coefficients)   DR7
;Band	    b	      [m(f/f0 = 0)]	m(f/f0 = 10b)
;u 	1.4 × 10-10	24.63	            22.12
;g 	0.9 × 10-10	25.11	            22.60
;r 	1.2 × 10-10	24.80	            22.29
;i 	1.8 × 10-10	24.36	            21.85
;z 	7.4 × 10-10	22.83	            20.32
;
;
;
;If clip set, negative fluxes are removed.
;The used index are returned 
;
;-


function asinh2pogs, magashin, band, clip=clip, used=used

;get the b param
case band of
    'u': b=1.4D-10
    'g': b=0.9D-10
    'r': b=1.2D-10
    'i': b=1.8D-10
    'z': b=7.4D-10
    else: return, -99
endcase

flux_f0=2*b*sinh(magashin/(-(2.5/alog(10.)))-alog(b))

if keyword_set(clip) then used=where(flux_f0 gt 0.) else $
  used=make_array(n_elements(flux_f0),/index)

pog=-2.5*alog10(flux_f0[used])

return, pog

end
