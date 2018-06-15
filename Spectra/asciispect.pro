;+
;
;  Take one of the final product of xidl spectra reduction 
;  and writes an ascii file that can be shared easily 
;  amoung non-xidl people.
;
;  fits     your fits file
;  outname  the filename to save
;  inflg    the type of spectrum to be read
;  fil_sig  error file
;
;  The output is meant to be the 'final' version of the spectrum.  
;  
;
;-


pro asciispect, fits, outname=outname, inflg=inflg, display=display, $
                fil_sig=fil_sig

;;set filename
if ~keyword_set(outname) then begin
   fits_pos=strpos(fits,'.fits')
   outname=strmid(fits,0,fits_pos)+'.dat'
endif


;;read the spectum with the appropriate inflg
flux=x_readspec(fits,INFLG=inflg,HEAD=head,SIG=sig,WAV=wav,fil_sig=fil_sig)

;;write the output out
forprint, wav, flux, sig, textout=outname, comm='Wave Flux Sigma'

;;if told, dispaly for comparison
if keyword_set(display) then begin
plot, wav, flux
oplot, wav, sig, color=250
x_specplot, fits, inflg=inflg
endif


end
