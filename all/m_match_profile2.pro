;+ 
; NAME:
; m_match_profile2 
;   Version 1.1
;
; PURPOSE:
;    Extract and calibrate HI and CO profiles using a same pattern of ellipses.
;    The CO center is selected by the user
; CALLING SEQUENCE:
;   m_profcalibH2, datainfo, astrohi, astroh2    
; INPUTS:
;   datainfo   ASCII file with information of the images to process: (HIimage, COimage, Xcenter (deg), Ycenter(deg),
;              STEP("),STOP("),PA(deg),HIPS(''/pix),COPS(''/pix),HIBUNIT, COBUNIT, HIBmin("),HIBmax("),COBmin("),COBmax("),HLUM(log Msun),
;              HIrms, COrms, FINAL, NSIGMA,INCL(deg))
;
;   astrohi    ASCII file with information of the HI astrometry (HICRVAL1,HICDELT1,HICRPIX1,HICRVAL2,HICDELT2,HICRPIX2)
;   astroco    ASCII file with information of the CO astrometry (COCRVAL1,COCDELT1,COCRPIX1,COCRVAL2,COCDELT2,COCRPIX2)
;
;
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


pro m_match_profile2,  datainfo, astrohi, astroh2
       if  N_params() LT 1  then begin 
          print, 'Syntax - ' +$
         'm_match_profile2,  datainfo, astrohi, astroh2 (v1.1)'
       return
       endif
 close, /all


;read data file
READCOL, datainfo, HIimage, COimage, Xcenter, Ycenter, STEPas, STOPas, PAdg, HIps,COps,HIBUNIT, COBUNIT, HIBmin, $
         HIBmax,COBmin,COBmax,HLUMlg,HIrms,COrms,FINALu,NSIG,INCLdg,FORMAT='A,A,F,F,F,F,F,F,F,A,A,F,F,F,F,F,F,F,A,F,F' ,/SILENT


READCOL, astrohi, NAME, HICRVAL1,HICDELT1,HICRPIX1,HICRVAL2,HICDELT2,HICRPIX2,FORMAT='A,D,D,D,D,D,D',/SILENT


READCOL, astroh2, NAME, COCRVAL1,COCDELT1,COCRPIX1,COCRVAL2,COCDELT2,COCRPIX2,FORMAT='A,D,D,D,D,D,D',/SILENT

close, /all

i=0

WHILE (i LT N_ELEMENTS(HIimage)) DO BEGIN


;read semimajor/minor axis

 fits = MRDFITS(COimage[i], 0, fitsheadHI)
  
  
  
  ;display image
  IF(N_ELEMENTS(fits[*,0]) LT 800) THEN RESCALE=0.5
  IF(N_ELEMENTS(fits[*,0]) GE 800 AND N_ELEMENTS(fits[*,0]) LT 1200) THEN RESCALE=1
  IF(N_ELEMENTS(fits[*,0]) GE 1200) THEN RESCALE=1.5
  
  Window,1,XSIZE=N_ELEMENTS(fits[*,0])/RESCALE,YSIZE=N_ELEMENTS(fits[0,*])/RESCALE
  disp=CONGRID(fits,N_ELEMENTS(fits[*,0])/RESCALE,N_ELEMENTS(fits[0,*])/RESCALE)
  TVSCL, disp 

  
  print, "Select galaxy center"
  CURSOR, Xc, Yc, /DEVICE, /DOWN 
  print, "Xc:",Xc*RESCALE," Yc:",Yc*RESCALE  
  print, "Select mjor axis direction"
  CURSOR, Xa, Ya, /DEVICE, /DOWN 
  print, "Xa:",Xa*RESCALE," Ya:",Ya*RESCALE  
  print, "Select minor axis direction"
  CURSOR, Xb, Yb, /DEVICE, /DOWN 
  print, "Xb:",Xb*RESCALE," Yb:",Ya*RESCALE  

  Apix=SQRT((Xc-Xa)^2+(Yc-Ya)^2)
  Bpix=SQRT((Xc-Xb)^2+(Yc-Yb)^2)
  
;extract HI profile

  Xcen=(Xc*RESCALE-COCRPIX1[i])*COCDELT1[i]+COCRVAL1[i]
  Ycen=(Yc*RESCALE-COCRPIX2[i])*COCDELT2[i]+COCRVAL2[i]

  
  HIxcen=(Xcen-HICRVAL1[i])/HICDELT1[i]+HICRPIX1[i]    
  HIycen=(Ycen-HICRVAL2[i])/HICDELT2[i]+HICRPIX2[i]    
  
  
   
  m_fixellipse, HIimage[i], Xcenter=HIxcen,Ycenter=HIycen,A=MAX([Apix,Bpix]),B=MIN([Apix,Bpix]),PA=PAdg[i], STEP=STEPas[i]/HIps[i],STOP=STOPas[i]/HIps[i]
   
  HIraw = strjoin([HIimage[i], '.raw'])
   m_profcalibHI, HIraw, PS=HIps[i], BUNIT=HIBUNIT[i], Bmin=HIBmin[i], Bmax=HIBmax[i], $
                  INCL=INCLdg[i], RMS=HIrms[i], NSIGMA=NSIG[i], FINAL=FINALu[i]

  
  
   
  
  COxcen=Xc*RESCALE   
  COycen=Yc*RESCALE   
   
   

  m_fixellipse, COimage[i], Xcenter=COxcen,Ycenter=COycen,A=MAX([Apix,Bpix]),B=MIN([Apix,Bpix]),PA=PAdg[i], STEP=STEPas[i]/COps[i],STOP=STOPas[i]/COps[i]
  COraw = strjoin([COimage[i], '.raw'])
  m_profcalibH2, COraw, PS=COps[i], BUNIT=COBUNIT[i], Bmin=COBmin[i], Bmax=COBmax[i], $
                 INCL=INCLdg[i], RMS=COrms[i], NSIGMA=NSIG[i], FINAL=FINALu[i], HLUM=HLUMlg[i]

  
  
  HIfinal = strjoin([HIraw, '.cal'])
  COfinal = strjoin([COraw, '.cal'])
   
  READCOL, HIfinal, HIRad, HIIoutput, HIIerr, SKIPLINE=0
  READCOL, COfinal, CORad, COIoutput, COIerr, SKIPLINE=0

  finalname = strjoin(['N',NAME[i], '_HICOcal.dat'])

  FORPRINT, HIRad, HIIoutput, HIIerr,COIoutput, COIerr, $
            COMMENT = "Rad(arcsec), HIImed , HIErrI, COImed , COErrI", TEXTOUT=finalname, /SILENT

 i=i+1
 
 PLOT, HIRad, HIIoutput, psym=4, yrange=[MIN([HIIoutput,COIoutput]),MAX([HIIoutput,COIoutput])]
 OPLOT,  HIRad, COIoutput, psym=2
 
 ENDWHILE

END
