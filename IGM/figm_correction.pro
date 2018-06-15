;Correct the flux  in one filter according to IGM blanket 
;This use the tabulated IGM correction from figm_calculator.

;lambda wavelenght array
;transm filter transmission
;zemi is the redshift 
;cIGM is output



PRO figm_correction, zemi, cIGM, FILT=filt

path=getenv("MIDL")+'/IGM/grid/'

;;open right grid
IF(filt EQ 'u_lris') THEN readcol, path+"IGM_U_LRIS_full.dat", redsh, IGM, skipline=1, /silent
IF(filt EQ 'b_lris') THEN readcol, path+"IGM_B_LRIS_full.dat", redsh, IGM, skipline=1, /silent
IF(filt EQ 'v_lris') THEN readcol, path+"IGM_V_LRIS_full.dat", redsh, IGM, skipline=1, /silent
IF(filt EQ 'r_lris') THEN readcol, path+"IGM_R_LRIS_full.dat", redsh, IGM, skipline=1, /silent
IF(filt EQ 'i_lris') THEN readcol, path+"IGM_I_LRIS_full.dat", redsh, IGM, skipline=1, /silent

IF(filt EQ 'v_lrisnew') THEN readcol, path+"IGM_V_LRISnew_full.dat", redsh, IGM, skipline=1, /silent
IF(filt EQ 'r_lrisnew') THEN readcol, path+"IGM_R_LRISnew_full.dat", redsh, IGM, skipline=1, /silent
IF(filt EQ 'i_lrisnew') THEN readcol, path+"IGM_I_LRISnew_full.dat", redsh, IGM, skipline=1, /silent


IF(filt EQ 'sd_lbc') THEN readcol, path+"IGM_SD_LBC_full.dat", redsh, IGM, skipline=1, /silent
IF(filt EQ 'b_lbc') THEN readcol, path+"IGM_B_LBC_full.dat", redsh, IGM, skipline=1, /silent
IF(filt EQ 'v_lbc') THEN readcol, path+"IGM_V_LBC_full.dat", redsh, IGM, skipline=1, /silent
IF(filt EQ 'r_lbc') THEN readcol, path+"IGM_R_LBC_full.dat", redsh, IGM, skipline=1, /silent
IF(filt EQ 'i_lbc') THEN readcol, path+"IGM_I_LBC_full.dat", redsh, IGM, skipline=1, /silent


IF(filt EQ 'f275w')  THEN readcol, path+"IGM_f275w_WFC3_full.dat", redsh, IGM, skipline=1, /silent
IF(filt EQ 'f336w')  THEN readcol, path+"IGM_f336w_WFC3_full.dat", redsh, IGM, skipline=1, /silent
IF(filt EQ 'f343n')  THEN readcol, path+"IGM_f343n_WFC3_full.dat", redsh, IGM, skipline=1, /silent
IF(filt EQ 'f390m')  THEN readcol, path+"IGM_f390m_WFC3_full.dat", redsh, IGM, skipline=1, /silent
IF(filt EQ 'f390w')  THEN readcol, path+"IGM_f390w_WFC3_full.dat", redsh, IGM, skipline=1, /silent
IF(filt EQ 'f438w')  THEN readcol, path+"IGM_f438w_WFC3_full.dat", redsh, IGM, skipline=1, /silent
  
;;get value
cIGM=interpol(IGM,redsh,zemi)


END
