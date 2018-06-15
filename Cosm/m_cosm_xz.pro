;+ 
; NAME:
; m_cosm_xz
;
; PURPOSE:
;    Calculate the cosmological distance X from z=0 to z=z
;
; CALLING SEQUENCE:
;   
;
; INPUTS:
;   z -- Redshift
;
; RETURNS:
;   x -- Cosmological pathlength
;
; OUTPUTS:
;
; OPTIONAL KEYWORDS:
;   H0 -- Hubbles constant (km/s/Mpc)
;   OM -- Omega Dark Matter
;   OV -- Lambda
;
; OPTIONAL OUTPUTS:
;
; COMMENTS:
;
; EXAMPLES:
;   
;
; PROCEDURES CALLED:
;  cosm_common
;  cosm_intxz
;  qromb
;
; REVISION HISTORY:
;   11-March-2004 Written by JXP
;-
;------------------------------------------------------------------------------

function m_cosm_intxz, z
common cosmolgy_cmmn, cosm_dm, cosm_k, cosm_h, cosm_Ob, cosm_L, cosm_r

  intxz = (1.+z)^2 / sqrt(cosm_L + (cosm_K)*(1+z)^2 + cosm_dm*(1+z)^3)
  return, intxz
end
  

function m_cosm_xz, z, H0=h0, OM=om, OV=ov, _EXTRA=extra

  common cosmolgy_cmmn

  if (N_params() LT 1) then begin 
    print,'Syntax - ' + $
             'xz = m_cosm_xz(z, H0=, OM=, OV=) [v1.1]'
    return, -1
  endif 

  m_cosm_common, H0=h0, Omegavac=OV, OmegaDM=OM, _EXTRA=extra

  xz = qromb('m_cosm_intxz', 0., z, /double)

  return, xz   ; Units are unknown

end

