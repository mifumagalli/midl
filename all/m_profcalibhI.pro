;+ 
; NAME:
; m_profcalibhI  
;   Version 1.1
;
; PURPOSE:
;    Calibrate HI surface brightness profiles at n*sigma
;
;
; CALLING SEQUENCE:
;   m_profcalibHI, filename, PS=, BUNIT=, Bmin=, Bmax=,  INCL=, RMS=, NSIGMA=, FINAL=
; INPUTS:
;   filename    ASCII file to process (output of m_fixellipse)
;   PS		Platescale
;   BUNIT       Raw unit ('JYMS','JYKS','KMS','KKMS')
;   Bmin        Beam FWHM minor in arcsec
;   Bmax        Beam FWHM major in arcsec  
;   RMS         Minimum sensitiviy (noise) 
;   FINAL       Desired final units ('MP' or 'N' for Msun/pc^2 or cm^-2)
;   NSIGMA	N*rms minimum sensitivity
;   INCL        Inclination angle (deg)  
;
; RETURNS:
;
; OUTPUTS:
; OPTIONAL OUTPUTS:
;
; COMMENTS:
;
; EXAMPLES:
;m_profcalibhI,  'N_6946_NA_mom0_THINGS.fits_prof_raw.dat', PS=1.5, BUNIT='JYMS', Bmin=5.61, Bmax=6.04, RMS=0.55, NSIGMA=10., INCL=40, FINAL='N'
;
; PROCEDURES/FUNCTIONS CALLED:
;  
;  
;
; REVISION HISTORY:
;   11-Nov-2008 Written by MF
;-
;------------------------------------------------------------------------------


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

pro m_profcalibhI, filename, PS=ps, BUNIT=bunit, Bmin=bmin, Bmax=bmax, NSIGMA=nsigma, INCL=incl, RMS=rms, FINAL=final
 
   if  N_params() LT 1  then begin 
         print, 'Syntax - ' +$
        'm_profcalibhI, filename, PS=, BUNIT=, Bmin=, Bmax=, NSIGMA=, INCL=, RMS=, FINAL= (v1.1)'
      return
   endif 

  close, /all

  readcol, filename, SMA, Intens, Npixel, SKIPLINE =1, /SILENT

  
  ;uncertainty calibration as in Wong&Blitz, 2002

   Pixelbeam=3.1416*(Bmin*Bmax)/(PS)^2
   SigmaIntens=RMS/SQRT(Npixel/Pixelbeam)
    
  ;calibration HI   (reference in Walter, 2008 THINGS; Fumagalli et al.,2008)
  
     IF (BUNIT EQ 'JYMS') THEN BEGIN
     Ndensity=1.1106E24*1E-3*(Intens)/(Bmin*Bmax)
     Nrms=1.1106E24*1E-3*(rms)/(Bmin*Bmax)
     Nerr=1.1106E24*1E-3*(SigmaIntens)/(Bmin*Bmax)
     ENDIF

     IF (BUNIT EQ 'JYKMS') THEN BEGIN
     Ndensity=1.1106E24*(Intens)/(Bmin*Bmax)
     Nrms=1.1106E24*(rms)/(Bmin*Bmax)
     Nerr=1.1106E24*(SigmaIntens)/(Bmin*Bmax)
     ENDIF

   
     IF (BUNIT EQ 'KMS') THEN BEGIN
     Ndensity=1.823E18*1E-3*Intens
     Nrms=1.823E18*1E-3*rms
     Nerr=1.823E18*1E-3*SigmaIntens
     ENDIF

     IF (BUNIT EQ 'KKMS') THEN BEGIN
     Ndensity=1.823E18*Intens
     Nrms=1.823E18*rms
     Nerr=1.823E18*SigmaIntens
     ENDIF

     IF (BUNIT NE 'JYMS' AND BUNIT NE 'JYKMS' AND BUNIT NE 'KMS' AND BUNIT NE 'KKMS') THEN BEGIN
     print, "I can't understand your units. Try with JYMS, JYKMS, KMS KKMS. No output!"
     ENDIF

  
  
   Mdensity=8.01E-21*Ndensity
   Mrms=8.01E-21*Nrms
   Merr=8.01E-21*Nerr
   Rad=SMA*PS
   
   IF (FINAL EQ 'N') THEN BEGIN
   Ioutput=Ndensity
   OutRMS=Nrms
   Ierr=Nerr
   ylab='HI (cm^-2)'
   ENDIF
      
   IF (FINAL EQ 'MP') THEN BEGIN
   Ioutput=Mdensity
   OutRMS=Mrms
   Ierr=Merr
   ylab='HI (Msun/pc^2)'
   ENDIF
   
   IF (FINAL NE 'MP' AND FINAL NE 'N') THEN BEGIN
   print, "No output units specified. Assuming column dnsity as default."
   Ioutput=Ndensity
   OutRMS=Nrms
   Ierr=Nerr
   ylab='HI (cm^-2)'
   ENDIF
  
  
   ;correction for inclination effect
   Ioutput=Ioutput*COS(INCL*3.1416/180)
   OutRMS=OutRMS*COS(INCL*3.1416/180)
   Ierr=Ierr*COS(INCL*3.1416/180)
   
   ;exclude value below N*sigma
   Ioutput=Ioutput[where(Ioutput GT NSIGMA*OutRMS)]
   Rad=Rad[where(Ioutput GT NSIGMA*OutRMS)]
   Ierr=Ierr[where(Ioutput GT NSIGMA*OutRMS)]
   
   ;plot, Rad, Ioutput, psym=6, xtitle="Rad (arcsec)", ytitle=ylab
   ;errplot, Rad, Ioutput-Ierr, Ioutput+Ierr
  
   file=strjoin([filename,'.cal'])
  
   print, "Results in ", file 
   FORPRINT, Rad, Ioutput, Ierr, COMMENT = "Rad(arcsec), Imed (ylab), ErrI" ,$
             TEXTOUT =file, /SILENT
  

end

