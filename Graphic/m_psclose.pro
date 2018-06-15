;+ 
; NAME:
; x_psclose
;   Version 1.1
;
; PURPOSE:
;    Call ps_close and reset a number of default values.  Point window
;   at X-windows.
;
; CALLING SEQUENCE:
;   
; INPUTS:
;
; RETURNS:
;
; OUTPUTS:
;
; OPTIONAL KEYWORDS:
;
; OPTIONAL OUTPUTS:
;
; COMMENTS:
;
; EXAMPLES:
;
; PROCEDURES/FUNCTIONS CALLED:
;  ps_close
;  set_plot
;
; REVISION HISTORY:
;   09-Feb-2004 Written by JXP
;-
;------------------------------------------------------------------------------

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

pro m_psclose, _extra=_extra, nohigh=nohigh

common rpsopen_files, filename1

sharpcorners, thick=!x.thick

;;Optional Keywords

  psclose
  set_plot, 'x'
  device, decompose=1
  !p.thick = 1
  !p.charthick = 1
  !p.font = -1
  !x.thick = 1
  !y.thick = 1

;;perform some operation on the output file 
ps_mod, filename1, high=(keyword_set(nohigh) eq 0), _extra=_extra
end
