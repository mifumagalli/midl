;+
; NAME:
;        NOBLANCK
;
; PURPOSE:
;
;        The NOBLANCK function returns a copy of Strings with all
;        blanck characters completely removed.
;
; CATEGORY:
;        String Processing.
;
; CALLING SEQUENCE:
;
;        Result = NOBLANCK( Strings )
;
; INPUTS:
;        Strings:  The string or array of strings to be processed.
;
; EXAMPLE:
;
; MODIFICATION HISTORY:
;        Written by:    MF, March 2009.
;-
function NOBLANCK, Str

;   Check integrity of input parameter

         NP        = N_PARAMS()
         if (NP ne 1) then message,'Must be called with 1 parameter, Str'

         sz        = SIZE(Str)
         ns        = n_elements(sz)
         if (sz(ns-2) ne 7) then message,'Parameter must be of string type.'
         ndim      = sz(0)

;   Find non-alphabetical characters

         strbyte   = BYTE(Str)
         return,strcompress(strbyte,/REMOVE_ALL)
end
