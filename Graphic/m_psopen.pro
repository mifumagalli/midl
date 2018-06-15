;+ 
; NAME:
; m_psopen
;   Version 1.1
;
; PURPOSE:
;  Calls ps_open and sets a number of default plotting values.  
;
; CALLING SEQUENCE:
;   
;   m_psopen, fil, /MAXS, /PORTRAIT
;
; INPUTS:
;   img       - Fits file or data
;
; RETURNS:
;   dat       - Data in fits file or the input data
;
; OUTPUTS:
;
; OPTIONAL KEYWORDS:
;  FSCALE      - Data is float
;
; OPTIONAL OUTPUTS:
;  HEAD        - Header
;
; COMMENTS:
;
; EXAMPLES:
;
;
; PROCEDURES/FUNCTIONS CALLED:
;
; REVISION HISTORY:
;   09-Feb-2004 Written by JXP. Modified by MF
;-
;------------------------------------------------------------------------------

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

pro m_psopen, psfile, MAXS=maxs, _EXTRA=extra

common rpsopen_files, filename1


if  N_params() LT 1  then begin 
    print,'Syntax - ' + $
      'm_psopen, psfile, /MAXS, _EXTRA= [V1.1]'
    return
endif 

if keyword_set(maxs) then begin
    YSIZ=7.9 ;8.25 
    XSIZ=10. ;11
endif


if keyword_set(color) then device, decompose=1

;; Optional Keywords
set_plot, 'x'
!p.font = 0
!p.charsize=1.8
device, decompose=0
if keyword_set(maxs) then psopen, PSFILE, /color, bpp=8, /inches, xsize=xsiz,$
  ysize=ysiz, _EXTRA=extra, /landsc else psopen, PSFILE, /color, bpp=8, _EXTRA=extra

!p.thick = 6
!x.thick = 8
!y.thick = 8
!p.charthick = 3
device, /times , /bold, /isolatin
;device, Set_Font='roman unicode', /TT_FONT

filename1=psfile

end
