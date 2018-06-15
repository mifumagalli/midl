;+ 
; NAME:
; m_plot_hico
;   Version 1.1
;
; PURPOSE:
; Plot CO and HI surface brightness profiles extracted with m_match_profile
; CALLING SEQUENCE:
;   m_plot_hico, name 
;
; INPUTS:
;  name   ASCII file with profile to plot
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
   
PRO m_plot_hico, name  

       if  N_params() LT 1  then begin 
          print, 'Syntax - ' +$
         'm_plot_hico, name (v1.1)'
       return
       endif

close, /all

  
  
  READCOL, name, HIRad, HIIoutput, HIIerr,COIoutput, COIerr

  
 PLOT, HIRad, HIIoutput, psym=4, yrange=[MIN([HIIoutput,COIoutput]),MAX([HIIoutput,COIoutput])]
 ERRPLOT, HIRad, HIIoutput-HIIerr, HIIoutput+HIIerr
 OPLOT,  HIRad, COIoutput, psym=2
 ERRPLOT, HIRad, COIoutput-COIerr, COIoutput+COIerr
 
end
