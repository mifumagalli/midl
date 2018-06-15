;+ 
; NAME:
; m_cosm_common
;
; PURPOSE:
;    Routine to initialize and set values in the Cosmology common
;    block named 'cosmolgy_cmmn'
;
; CALLING SEQUENCE:
;   cosm_common
;
; INPUTS:
;   H0 =   Hubbles constant in km/s/Mpc
;   Omegavac = Lambda value   [Default: 0.7]
;   OmegaDM = Omega for Dark Matter  [Default: 0.3]
;
; RETURNS:
;
; OUTPUTS:
;   
; OPTIONAL KEYWORDS:
;  /W06MAP -- Applies the WMAP cosmolgy from 2006
;  /SILENT -- No screen output
;
; OPTIONAL OUTPUTS:
;
; COMMENTS:
;
; EXAMPLES:
;
; PROCEDURES CALLED:
;
; REVISION HISTORY:
;   22-Nov-2003 Written by JXP
;-              Inlcuded WMAP5,BOLSHOI by MF
;------------------------------------------------------------------------------

pro m_cosm_common, H0=h0, Omegavac=omegavac, OmegaDM=omegaDM, WMAP7=WMAP7, WMAP5=WMAP5,$
                 WMAP3=WMAP3, WMAP9=WMAP9, P13=P13, BOLSH=BOLSH, DEF=def, SILENT=silent, TEST=TEST

common cosmolgy_cmmn, cosm_dm, cosm_K, cosm_h, cosm_Ob, cosm_L, cosm_r, sigma_8
common cosmolgy_other, omega_baryon

;common cosmolgy_cmmn, cosm_dm, cosm_K, cosm_h, cosm_L, cosm_r


  ;;This set the default cosmology (72,0.7,0.3)
  if not keyword_set( OmegaDM ) then cosm_dm = 0.3 else cosm_dm = omegaDM
  if size( H0 , /type) EQ 0 then cosm_h = 72. else cosm_h = H0
  if size( Omegavac, /type) EQ 0 then cosm_L = 0.7 else cosm_L = Omegavac
  sigma_8=0.8
  omega_baryon=0.0455

  ;;redefine if told to do so
  if keyword_set(P13) then begin
      cosm_h = 67.8
      cosm_dm = 0.308
      cosm_L = 1-cosm_dm
      sigma_8= 0.823
      omega_baryon= 0.0461 ;;Need to check this one
      if not keyword_set(SILENT) then print, "Using Planck13"
  endif

  ;;redefine if told to do so
  if keyword_set(WMAP9) then begin
      cosm_h = 69.7
      cosm_dm = 0.236
      cosm_L = 0.764
      sigma_8= 0.817
      omega_baryon= 0.0461
      if not keyword_set(SILENT) then print, "Using WMAP9"
  endif
 
  ;;corrected to recommanded values Sep 2012
  if keyword_set(WMAP7) then begin
      cosm_h = 70.4
      cosm_dm = 0.272
      cosm_L = 0.728
      sigma_8=0.809
      omega_baryon=0.0456
      if not keyword_set(SILENT) then print, "Using WMAP7"
  endif
  
  if keyword_set(WMAP5) then begin
      cosm_h =70.5
      cosm_dm =0.27
      cosm_L =0.73      
      sigma_8=0.81      
      if not keyword_set(SILENT) then print, "Using WMAP5"
  endif

  if keyword_set(WMAP3) then begin
      cosm_h = 73.
      cosm_dm = 0.24
      cosm_L = 0.76
      sigma_8=0.80      
      if not keyword_set(SILENT) then print, "Using WMAP3"
  endif

  if keyword_set(BOLSH) then begin
      cosm_h = 70.
      cosm_dm = 0.27
      cosm_L = 0.73
      sigma_8=0.82      
      if not keyword_set(SILENT) then print, "Using BOLSHOI values"
  endif

  if keyword_set(TEST) then begin
      cosm_h = 70.
      cosm_dm = 0.3
      cosm_L = 0.7
      sigma_8=0.
      if not keyword_set(SILENT) then print, "Using TEST"
  endif

  if keyword_set(DEF) then begin
      cosm_h = 72.
      cosm_dm = 0.3
      cosm_L = 0.7
      sigma_8=0.80
      if not keyword_set(SILENT) then print, "Using (72,0.7,0.3)"
  endif




  cosm_K = 1.d - cosm_dm - cosm_L
  cosm_r = 2.4725753d-05
  cosm_Ob = 0.
  if not keyword_set(SILENT) then begin
      print, 'm_cosm_common: Using this cosmology --'
      print, 'Omega_m = ', cosm_dm
      print, 'H (km/s/Mpc) = ', cosm_h
      print, 'Omega_L  = ', cosm_L
      print, 'sigma_8  = ', sigma_8
      print, 'Omega_bar  = ', omega_baryon
  endif     

  return

end

