;+
;
;   Compute the Mstar mag following equation 6 of Sbarufatti 2006
; 
;  redshift  the redshift array 
;  Mrstar    the host galaxy mag at z=0 
;  mstarz    Mstar in output
;
;-


pro redlim_mstar, redshift, Mrstar, mstarz=mstarz

  ;;correct for evolution with Bressan 1998 model (1e12)
  readcol, getenv("MIDL")+"/Spectra/blazar/bressan1998.txt", mod_time,$
           mod_mv, mod_vr, /sil
  mod_mr=mod_mv-mod_vr 
  mod_mr=mod_mr-max(mod_mr)+Mrstar  ;;normalize to given mag

  ;;correct time for new cosmology 
  mod_time=mod_time+13.743567-max(mod_time)
  mod_z=m_cosm_ztime((13.743567-mod_time)*1d9>0.,/wmap7)
  
  ;;find evolution corrected magnitude
  mstarz=interpol(mod_mr,mod_z,redshift)
  
end
