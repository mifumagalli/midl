;+
;PURPOSE
;	convert a spectrum to an sdss magnitude... note always does
;	AB magnitudes... none of that asinh stuff
;SYNTAX
;	mags=spec2mag(wave, flux[, totflux=totflux, filtfile=filtfile,
;			filtwave=filtwave, filt_thru=filt_thru, spec=spec])
;INPUTS
;	wave: wavelength in angstroms
;	flux: flux in 10^(-17) ergs/s/cm^2/ang
;	filtfile: if you don't want sdss magnitude then provide a filter file
;		of the filter you want
;	filtwave: filter wavelength (alternative to filtfile inputs)
;	filt_thru: filter thru-put (alternative to filtfile inputs)
;       spect: A fits file of a spectrum supported by x_readspec.
;             (alternative to wave flux inputs)  
;       inflg: keyword associated to spect for x_readspec
;
;KEYWORD 
;        displ  if set, a plot is produced       
;
;OUTPUTS
;	mags: 5 sdss magnitudes
;	totflux: the total flux prior to conversion to magnitudes
;NOTES: 
;	portions hacked from filter_thru.pro
;	only calculates the pogson magnitude. this version doesn't deal with 
;	that asinh mess
;Written by R. da Silva, UCSC
;Optimized for costum filters and spectra by MF, Apr 2010
;-

function spec2mag, wave, flux, totflux=totflux, filtfile=filtfile, $
	filtwave=filtwave, filt_thru=filt_thru, spect=spect, inflg=inflg, displ=displ


if keyword_set(spect) then flux=x_readspec(spect,wav=wave,inflg=inflg)


c=double(2.9979000e+10)
wave=double(wave)
flux=double(flux)

;set default sdss file
if (not keyword_set(filtwave) or not keyword_set(filtfile)) then begin
    filter_prefix='sdss_jun2001'
    ffiles = [filter_prefix+'_u_atm.dat', filter_prefix+'_g_atm.dat', $
              filter_prefix+'_r_atm.dat', filter_prefix+'_i_atm.dat', $
              filter_prefix+'_z_atm.dat']
endif 

if keyword_set(filtfile) then ffiles=filtfile

if keyword_set(filtwave) then ffiles=1



nfiles=n_elements(ffiles)
;lambda=float([3543,4770,6231,7625,9134])   
mags=dblarr(nfiles)
totflux=dblarr(nfiles)

for ifile=0, nfiles-1 do begin
   if not keyword_set(filtfile) then begin
      filename = filepath(ffiles[ifile], $
                          root_dir=getenv('IDLUTILS_DIR'), $
                          subdirectory=['data','filters'])
   endif else filename=filtfile[ifile]
   if not keyword_set(filtwave) then begin
      readcol, filename, fwave, fthru, /silent 
   endif else begin
      fwave=filtwave
      fthru=filt_thru
   endelse
                                ;get central wavelenght of the filter
    lambda=tsum(fwave, fthru*fwave)/tsum(fwave, fthru)
    
    linterp, fwave, fthru, wave, interpfilt, missing=0
                                ;added this line aboue missing
                                ;had to convert everything to doubles
    totfilt=tsum(wave, interpfilt)
    totflux[ifile]=1d-17*lambda^2*1d-8*$
      tsum(wave, flux*interpfilt)/totfilt/c
        ; flux convolved with filter/ integral over filter
        ;times lambda^2/c
        ;factor of 1d-17 for the initial units
	;factor of 1d-8 to convert from angstroms

    if keyword_set(displ) then begin
        !p.multi=[0,2,2]
        plot, wave, flux, title='Input spectrum'
        plot, fwave, fthru, title='Filter'
        plot, wave, flux*interpfilt/totfilt,  title='Filter*Spectrum'
        !p.multi=0
    endif
    
endfor

mags=-2.5*alog10(totflux)-48.6

;if total(finite(mags)) ne 5 then begin
;    splog, 'Mag calculated to not be real'
;    stop
;    return, -1
;endif

return, mags
end
