;This procedure compute the SFR for an objetc with known flux  

;flux          --> the object flux in cgs (or mag in AB mag) 
;error         --> the error associated to the flux (or mag in AB mag) 
;redshift      --> object redshift
;TYPE          --> The calibration. Currently implemented MAD1500
;IGM           --> If set to a filter known to figm_correction, compute the IGM correction
;dust          --> If set, a dust correction is computed



PRO compute_sfr, flux_i, error_i, redshift, type=type, igm=igm, dust=dust


;check if it is a mag
if (flux_i gt 1) then begin

    splog, 'Convert mag to flux'
;find flux under AB mag assumption
    flux=10^(-0.4*(flux_i+48.6))
    SN=1.0857/error_i
    error=flux/SN
endif else begin
    flux_i=flux
    error_i=error
endelse



;compute distance
distance_calculator, redshift, Distance, /lum, /DEF
cgs, cns
Distance=Distance*1D6*cns.pc


;Luminosity + simple K-correction
Lum=(flux*4*!PI*Distance^2)/(1+redshift)
ErrLum=(error*4*!PI*Distance^2)/(1+redshift)


;Star formation 
;MAD1500: Madau (1998) at 1500A, corrected for Chabrier IMF (/1.59) from Salim, 2007


case type of 
    'MAD1500' : sfrcalib=7.9114D-29
    else : splog, "SFR unknown!"
endcase
splog, "SFR calibration ", sfrcalib
SFR=Lum*sfrcalib
errSFR=ErrLum*sfrcalib
splog, 'Calibrated SFR (Uncorrected): ', SFR
splog, 'Error Calibrated SFR (Uncorrected): ', errSFR

if keyword_set(IGM) then begin
;IGM correction 
    figm_correction, redshift, cIGM, FILT=IGM
    splog, "IGM correction: ", cIGM
endif else cIGM=1.


if keyword_set(dust) then begin

;dust reddy 2004
    if(redshift lt 2) then fdust=4.3
    if(redshift ge 2 and redshift lt 2.5) then fdust=4.4
    if(redshift ge 2.5) then fdust=4.7
    
    cdust=fltarr(n_elements(sfr))
    
    highsf=where(sfr*cigm gt 20, nhsf, complement=lowsfr)
    if(nhsf gt 0) then cdust[highsf]=fdust
    if(nhsf lt n_elements(sfr)) then cdust[lowsfr]=2.3

    splog, 'Dust correction: ', cdust
    
endif else cdust=1.

;sfr corrected
sfrcorr=sfr*cdust*cIGM
SFRerrcorr=errSFR*cdust*cIGM


splog, 'Corrected SFR: ', sfrcorr
splog, 'Error Corrected SFR: ', SFRerrcorr




END
