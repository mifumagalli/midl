;This function return the EBV dust from Schlegel map.
;If filter keyword is set, it computes the extinction according to 
;standard extinction law.

;RA      --> input RA hh:mm:ss
;DEC     --> input DEC dd:mm:ss
;A       --> if band keyword set, it returns the extinction 
;BAND    --> set to a filter name to get the extinction back
;GG      --> set to pass the rag and deg instead of RA-DEC in hh
;LAMBDA  --> if set to a wavelenght in A, uses MW extinction 
;            law direclty (band is overwritten)   
;FILTER  --> if filter id provided, return Ax with convolution
;            rescaling as in Schlafly & Finkbeiner 2011

function ebv_dust, ra, dec, a=a, band=band, gg=gg, lambda=lambda, filter=filter


  ;;convert in deg-rag - galactic
  if ~keyword_set(gg) then x_radec, ra, dec, rag, deg $
  else begin
     rag=ra
     deg=dec
  endelse
  euler, rag, deg, l, b, 1
  

  
  EBVdust=dust_getval(l,b,ipath=getenv('MIDL')+$
                      '/Imaging/Photometry/dustmap/',/interp) 

  ;;Use tabulated values
  if keyword_set(band) then begin
     ;;Optical from Cardelli, 1989
     ;;FUV/NUV galex from Wyder et al. 2005, ApJ, 619, L15
     ;;A(Ha)=0.6*Ab from Kennicutt, 2008    
     ;;SDSS is computed using central wavelenght and MW extinction law
     case band of
        'FUV': A=8.376*EBVdust
        'NUV': A=8.741*EBVdust
        'U': A=4.8*EBVdust
        'B': A=4.1*EBVdust
        'V': A=3.1*EBVdust
        'R': A=2.3*EBVdust
        'I': A=1.5*EBVdust
        'HA': A=0.6*4.1*EBVdust
        'usdss': A=4.9*EBVdust
        'gsdss': A=3.7*EBVdust
        'rsdss': A=2.7*EBVdust
        'isdss': A=2.0*EBVdust
        'zsdss': A=1.45*EBVdust
        else: print, 'Filter not found' 
     endcase 
     
  endif

  
  ;;If frequency specified, uses MW extinction law
  if keyword_set(lambda) then begin
     A_v=3.1*EBVdust
     A=dust_law(lambda,A_v)
  endif


  ;;if filter specified
  if keyword_set(filter) then begin
     
     ;;make wavearray
     wavex=mkarr(1500.,15000.,1.)
     nwavx=15000-6500 ;index at 1mu
    
     ;;get o'donnell
     odonAlam=ext_odonnell(wavex,3.1)
     ;;get fitzpatrick
     fitzAlam=x_alav(wavex,RV=3.1)
     
     ;;renorm to 1micron
     fitzAlam=fitzAlam/fitzAlam[nwavx]
     
     ;;get Sch at 1 mu with rescale 
     EBVdust=EBVdust*0.78
     A_v_odo=3.1*EBVdust
     a1mu=odonAlam[nwavx]*A_v_odo
     A_l_fitz=a1mu*fitzAlam
     

     ;;load filter 
     getfilter, filter, lambda=flambda, trans=ftrans
     flambda=flambda[where(ftrans gt 0.001)]
     ftrans=ftrans[where(ftrans gt 0.001)]


     alfitint=interpol(A_l_fitz,wavex,flambda)
     deltal=flambda-shift(flambda,1)
     deltal[0]=0

     ;;make the integral
     A=-2.5*alog10(total(ftrans*wavex^(-2)*10^(-alfitint/2.5)*deltal)/total(ftrans*wavex^(-2)*deltal))
     
     plot, wavex, A_l_fitz, xrange=[1500,15000]
     oplot, flambda, ftrans/max(ftrans)*max(A_l_fitz), line=1
     oplot, flambda, alfitint, color=250
     
  endif
  
return, EBVdust

end
