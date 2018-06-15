;+ 
; NAME:
; m_stdmag   
;   Version 1.1
;
; PURPOSE:
;    Does aperture photometry on the objects in the std_* lists.
;   Drives the program x_aper
;
; CALLING SEQUENCE:
;   m_stdmag, struct, SKYR=, STRROOT=, OUTROOT=
;
; INPUTS:
;   struct -- dimg_strct defining the images of interest.  This
;             program focuses on the STD frames only.
;
; RETURNS:
;
; OUTPUTS:
;
; OPTIONAL KEYWORDS:
;  SKYR =  2 float array describing sky annulus in arcsec
;            (default=15,30)
;  STRROOT=  Path for standard stars [default: 'Photo/std_']
;  OUTROOT=  Path for output of magnitudes [default: 'Photo/mag_']
;
; OPTIONAL OUTPUTS:
;
; COMMENTS:
;
; EXAMPLES:
;   xdimg_stdmag, dimg
;
;
; PROCEDURES/FUNCTIONS CALLED:
;  X_APER
;  X_CCDINF
;
; REVISION HISTORY:
;   07-Aug-2001 Written by JXP. Modified by MF
;-
;------------------------------------------------------------------------------


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

pro m_stdmag, struct, SKYR=skyr, STRROOT=strroot, OUTROOT=outroot

;
  if  N_params() LT 1  then begin 
      print, 'Syntax - ' +$
        'm_stdmag, struct, SKYR=, STRROOT=, OUTROOT= (v1.1)'
      return
  endif 

  close, /all

;  Optional Keywords

  if not keyword_set( SKYR ) then SKYR = [15., 30.]
  if not keyword_set( STRROOT ) then strroot='Photo/std_'
  if not keyword_set( OUTROOT ) then outroot='Photo/mag_'


;  Find the Standard Stars

  stds = where(struct.type EQ 'STD' AND struct.flg_anly NE 0, nstds)
  if nstds EQ 0 then begin
      print, 'No standard stars to analyse!'
      return
  endif

;  Set aperture

  
  apr = 7.0
  skyrad = [12.0, 30.0]

; Set files
  datfil = ' '
  listfil = ' '
  nm = ''

;  Loop on all frames
   
    
 FOR q=0,nstds-1 do begin

      ; Fixed a bug.
      ; The std .dat file is really dealing with the FINAL image 
      ; and so it needs to look for the f_ prefix.
      strfil = strjoin([strroot,strmid(strtrim(struct[stds[q]].img_root,2),0, $
                              strlen(strtrim(struct[stds[q]].img_root,2))-5), '.dat'])
       
      
      ; a = findfile(strfil, count=cnt) findfile is obsolete.
       a = file_search(strfil, count=cnt)
      if cnt EQ 0 then continue

      ; Data

      data = xmrdfits(strtrim(struct[stds[q]].img_final,2), 0, head, /fscale, /silent)
    
      ; Bad pixels (saturation)
      badpix = [-30000., struct[stds[q]].satur]

      nlin = numlines(strfil)
      nam = strarr(nlin-2)
      xc = fltarr(nlin-2)
      yc = fltarr(nlin-2)

      openr, 1, strfil
      readf, 1, datfil
      readf, 1, listfil
      for j=0,nlin-3 do begin
          readf, 1, FORMAT='(a15,f,f)', nm, x, y
          nam[j] = nm
          xc[j] = x
          yc[j] = y
     endfor
   

      close, 1

      ; Aperture Photometry
      ; Be careful here!
      ; Give x_aper a gain of 1 if we already applied gain in xdimg_proc
      ; Give x_aper a gain of struct[stds[q]].gain if we didn't apply gain.
     
      ;set platescale
      IF (struct[stds[q]].filter EQ 'U') THEN BEGIN
      pixapr=apr/0.135
      pixskyrad = skyrad/0.135
      ENDIF
      
      IF (struct[stds[q]].filter NE 'U') THEN BEGIN
      pixapr=apr/0.211
      pixskyrad = skyrad/0.211
      ENDIF
       
      x_aper, data, xc, yc, mags, errap, sky, skyerr, struct[stds[q]].gain, $
        pixapr, pixskyrad, badpix, /silent, /exact
      
         
      ; Deal with exposure time [After to get correct errors]
          
      for j=0,nlin-3 do mags[j] = mags[j] + 2.5*alog10(struct[stds[q]].exp)
  
      ; Output
      outfil = strjoin([outroot,'f_',strmid(strtrim(struct[stds[q]].img_root,2),0, $
                                       strlen(strtrim(struct[stds[q]].img_root,2))-5), '.dat'])
      openw, 1, outfil
      printf, 1, struct[stds[q]].filter
      printf, 1, listfil
      printf, 1, struct[stds[q]].AM

      for j=0,nlin-3 do printf, 1, FORMAT='(a15,f,f)', nam[j], mags[j], errap[j]
      close, 1
      
      ; Release memory
      delvarx, nam, xc, yc, sky, skyerr, mags, errap
  endfor

  print, 'All done with photometry!'
  
end
