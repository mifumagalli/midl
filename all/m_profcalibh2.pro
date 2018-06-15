;+ 
; NAME:
; m_profcalibh2  
;   Version 1.1
;
; PURPOSE:
;    Calibrate H2 surface brightness profiles at n*sigma

;
; CALLING SEQUENCE:
;   m_profcalibH2, filename, PS=, BUNIT=, Bmin=, Bmax=, HLUM=, INCL=, RMS=, FINAL=, NSIGMA=
; INPUTS:
;   filename    ASCII file to process (output of m_fixellipse)
;   PS		Platescale
;   BUNIT       Raw unit ('JYMS','JYKS','KMS','KKMS')
;   Bmin        Beam FWHM minor in arcsec
;   Bmax        Beam FWHM major in arcsec  
;   RMS         Minimum sensitiviy (noise) 
;   FINAL       Desired final units ('MP' or 'N' for Msun/pc^2 or cm^-2)
;   HLUM	H band luminosity (Log)
;   INCL        Inclination angle (deg)  
;   NSIGMA      Number of sigma
; RETURNS:
;
; OUTPUTS:
;  
;
; OPTIONAL OUTPUTS:
;
; COMMENTS:
;
; EXAMPLES:
;m_profcalibH2,  filename, PS=1.5, BUNIT='JYMS', HLUM=10.25, Bmin=5.61, Bmax=6.04, RMS=0.55, INCL=40, FINAL='N'
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

pro m_profcalibh2, filename, PS=ps, BUNIT=bunit, Bmin=bmin, Bmax=bmax, HLUM=hlum, INCL=incl, RMS=rms, FINAL=final, NSIGMA=nsigma
 
   if  N_params() LT 1  then begin 
         print, 'Syntax - ' +$
        'm_profcalibH2, filename, PS=, BUNIT=, Bmin=, Bmax=, HLUM=, INCL=, RMS=, FINAL=, NSIGMA= (v1.1)'
      return
   endif 

  close, /all

  readcol, filename, SMA, Intens, Npixel, SKIPLINE =1, /SILENT

   ;uncertainty calibration as in Wong&Blitz, 2002

   Pixelbeam=3.1416*(Bmin*Bmax)/(PS)^2
   SigmaIntens=double(RMS/SQRT(Npixel/Pixelbeam))
   
   ; X conversion factor Boselli, 2002
   
   XLOG=-0.38*HLUM+24.23

  ;calibration CO
       
     IF (BUNIT EQ 'JYMS') THEN BEGIN
     Ndensity=10^XLOG*1E-3*Intens*92.5486/(Bmin*Bmax)
     Nrms=10^XLOG*1E-3*rms*92.5486/(Bmin*Bmax)
     Nerr=10^XLOG*1E-3*SigmaIntens*92.5486/(Bmin*Bmax)
     ENDIF
     
     IF (BUNIT EQ 'JYKMS') THEN BEGIN
     Ndensity=10^XLOG*Intens*92.5486/(Bmin*Bmax)
     Nrms=10^XLOG*rms*92.5486/(Bmin*Bmax)
     Nerr=10^XLOG*SigmaIntens*92.5486/(Bmin*Bmax)
     ENDIF

     IF (BUNIT EQ 'KMS') THEN BEGIN
     Ndensity=10^XLOG*1E-3*Intens
     Nrms=10^XLOG*1E-3*rms
     Nerr=10^XLOG*1E-3*SigmaIntens
     ENDIF

     IF (BUNIT EQ 'KKMS') THEN BEGIN
     Ndensity=10^XLOG*Intens
     Nrms=10^XLOG*rms
     Nerr=10^XLOG*SigmaIntens
     ENDIF

     IF (BUNIT NE 'JYMS' AND BUNIT NE 'JYKMS' AND BUNIT NE 'KMS' AND BUNIT NE 'KKMS') THEN BEGIN
     print, "I can't understand your units. Try with JYMS, JYKMS, KMS KKMS. No output!"
     ENDIF

  
  
   Mdensity=1.602E-20*Ndensity
   Mrms=1.602E-20*Nrms
   Merr=1.602E-20*Nerr
   Rad=SMA*PS
   
   IF (FINAL EQ 'N') THEN BEGIN
   Ioutput=Ndensity
   OutRMS=Nrms
   Ierr=Nerr
   ylab='H2 (cm^-2)'
   ENDIF
      
   IF (FINAL EQ 'MP') THEN BEGIN
   Ioutput=Mdensity
   OutRMS=Mrms
   Ierr=Merr
   ylab='H2 (Msun/pc^2)'
   ENDIF
   
   IF (FINAL NE 'MP' AND FINAL NE 'N') THEN BEGIN
   print, "No output units specified. Assuming column density as default."
   Ioutput=Ndensity
   OutRMS=Nrms
   Ierr=Nerr
   ylab='H2 (cm^-2)'
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

