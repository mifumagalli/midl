;+
;
;
;  Take as an input the filter transmission curves in lambda space
;  with wavelengh in AA, the CCD quantum efficiency and the atmospheric extinction and generates
;  a total transmission filter. All the units in lambda should be in AA.
;  The code gen_trans.py  renormalize the integral, so it is not a problem to have
;  relative units in the transmission.
;
;  filt_tra   --> name of the filter transmission
;  atm_ext    --> name of the amospheric transmission
;  ccd_qe     --> name of the ccd QE
;  path       --> where the above files are located  
;  name       --> name to assign to the filter
;
;-




pro total_transmission, filt_tra, atm_ext, ccd_qe, path=path, name=name

if not keyword_set(path) then path='./'
if not keyword_set(name) then path='tottrans'


;;read the data
readcol, path+filt_tra, l_f, t_f, /sil, format='D,D'
readcol, path+atm_ext, l_a, t_a, /sil, format='D,D'
readcol, path+ccd_qe, l_q, t_q, /sil, format='D,D'

;;rewrite a tmp file
forprint, l_f, t_f,  textout='tmp_bpz_'+filt_tra, /nocom
forprint, l_a, t_a,  textout='tmp_bpz_'+atm_ext, /nocom
forprint, l_q, t_q,  textout='tmp_bpz_'+ccd_qe, /nocom

;;run filter convolution
spawn, 'python '+getenv("BPZPATH")+'/gen_trans.py '+name+' ./tmp_bpz_'+filt_tra+' ./tmp_bpz_'+atm_ext+$
  ' ./tmp_bpz_'+ccd_qe

;;read new filter
readcol, name+'.res', new_fl, new_ft, /sil


;;some display
!P.MULTI=[0,2,2]

plot, l_a, t_a, xtitle='Lambda', ytitle='Atmosphere'

plot, l_q, t_q, xtitle='Lambda', ytitle='CCD QE'

plot, l_f, t_f, xtitle='Lambda', ytitle='Filter transmission', yrange=[0,1]
oplot, new_fl, new_ft, line=1

!P.MULTI=[0,0,0]


;;clean 
spawn, 'rm -f tmp_bpz_*'
spawn, 'rm -f plot_trans.sm'

end
