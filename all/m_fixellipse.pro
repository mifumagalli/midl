;+ 
; NAME:
; m_fixellipse  
;   Version 1.1
;
; PURPOSE:
;    Extract surface brightness profiles in elliptical annuli
;
;
; CALLING SEQUENCE:
;   m_fixellipse, filename, Xcenter=, Ycenter=, A=, B=, PA=, STEP=, STOP=
; INPUTS:
;   filename    fits file to process
;   Xcenter     the X center position in pixel
;   Ycenter     the Y center position in pixel 
;   A           the semimajor axis in pixel 
;   B           the semiminor axis in pixel
;   PA          PA in degree (+ counterclockwise)    
;   STEP        the increment along A in pixel
;   STOP        the maximum radius at which is stops
;
; RETURNS:
;
; OUTPUTS:
; A file filename_profile which contains for each position the mean intensity.
; OPTIONAL KEYWORDS:
;
; OPTIONAL OUTPUTS:
;
; COMMENTS:
;
; EXAMPLES:
;m_fixellipse, "N_6946_NA_mom0_THINGS.fits",Xcenter=500,Ycenter=500,A=400,B=400,PA=1, STEP=60,STOP=300 
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

pro m_fixellipse, filename, Xcenter=XC, Ycenter=YC, A=a, B=b, PA=pa, STEP=step, STOP=stop

   if  N_params() LT 1  then begin 
         print, 'Syntax - ' +$
        'm_fixellipse, filename, Xcenter=, Ycenter=, A=, B=, PA=, STEP=, STOP= (v1.1)'
      return
   endif 

  close, /all

  fits = MRDFITS(filename, 0, fitshead, /dscale, /silent)
  
  
  
  ;display image
  IF(N_ELEMENTS(fits[*,0]) LT 800) THEN RESCALE=0.5
  IF(N_ELEMENTS(fits[*,0]) GE 800 AND N_ELEMENTS(fits[*,0]) LT 1500) THEN RESCALE=1
  IF(N_ELEMENTS(fits[*,0]) GE 1500) THEN RESCALE=1.5
  
  TVLCT, 0, 0, 0, 100
  Window,1,XSIZE=N_ELEMENTS(fits[*,0])/RESCALE,YSIZE=N_ELEMENTS(fits[0,*])/RESCALE
  disp=CONGRID(fits,N_ELEMENTS(fits[*,0])/RESCALE,N_ELEMENTS(fits[0,*])/RESCALE)
  TVSCL, disp 
  TVLCT, 255, 255, 0, 100 
  
  
   ;define ellipse out 
   i=0
   e=SIN(ACOS(1.*B/A)) 
   PArad=PA*3.1416/180
   
   Imean=fltarr(STOP/(1.*STEP)+1)
   SMA=fltarr(STOP/(1.*STEP)+1)
   npixel=fltarr(STOP/(1.*STEP)+1)
   
   
   WHILE (STEP*i LT STOP) DO BEGIN 
     
   Aout=STEP*i+STEP
   Bout=Aout*SQRT(1-e^2)
   F1out=Aout*e
   F2out=-F1out
   XF1out=F1out*COS(PArad)
   XF2out=F1out*SIN(PArad)
   YF1out=F2out*COS(PArad)
   YF2out=F2out*SIN(PArad)
   
   Ain=STEP*i 
   Bin=Ain*SQRT(1-e^2)
   F1in=Ain*e
   F2in=-F1in
   XF1in=F1in*COS(PArad)
   XF2in=F1in*SIN(PArad)
   YF1in=F2in*COS(PArad)
   YF2in=F2in*SIN(PArad)
   
   ;display ellipse
   TVELLIPSE, Ain/RESCALE, Bin/RESCALE, XC/RESCALE, YC/RESCALE,  PA, COLOR=100
   
   ;find points within annulus    
   Xmax=N_ELEMENTS(fits[*,0])
   Ymax=N_ELEMENTS(fits[0,*])
   Index=DINDGEN(Xmax,Ymax)
   
   
    DISTFOCin=sqrt((FLOOR(Index/Xmax)-XC-XF1in)^2+((Index-Xmax*FLOOR(Index/Xmax))-YC-YF1in)^2)+ $
             sqrt((FLOOR(Index/Xmax)-XC-XF2in)^2+((Index-Xmax*FLOOR(Index/Xmax))-YC-YF2in)^2)
    DISTFOCout=sqrt((FLOOR(Index/Xmax)-XC-XF1out)^2+((Index-Xmax*FLOOR(Index/Xmax))-YC-YF1out)^2)+ $
              sqrt((FLOOR(Index/Xmax)-XC-XF2out)^2+((Index-Xmax*FLOOR(Index/Xmax))-YC-YF2out)^2)
   
   ;compute mean intensity in annuli masking negative pixel
   index = WHERE(DISTFOCin GT 2*Ain AND DISTFOCout LE 2*Aout AND fits GE 0., count)
   IF count NE 0 THEN npixel[i]=N_ELEMENTS(fits[index]) ELSE npixel[i]=0.
   IF count NE 0 THEN Imean[i]=MEAN(fits[index],/DOUBLE) ELSE Imean[i]=0.
   SMA[i]=Aout   
   
   i=i+1  
     
   ENDWHILE   	
   
   close, /all
   
   ;write data to file
   outfil = strjoin([filename, '.raw'])
   FORPRINT, SMA[where(SMA NE 0.)], Imean[where(SMA NE 0.)], npixel[where(SMA NE 0.)], COMMENT = "SMA, Imed, Npixel" ,$
             TEXTOUT =outfil, /SILENT
    

   print, 'Profile extraction ', filename, ' done'


end
 
