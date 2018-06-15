;procedure that plot the sed/photo-z results of one object as computed by eazy

;objid       --->    id of the object to fit
;[filters]   --->    (optional) array with code of filters used in the fit, 
;                   as in master.FILTER.RES.info. To overplot filters.
;[path]      --->    (optional) path to access files 
;[print]     --->    set keyword to save file
;[savename]  --->    name to eps file to save (no extension!)
;[yrange]    --->    if set, change yrange

;e.g. for LRIS_old [U,V,R,I]--->[127,128,129,130]



PRO eazy_plot, objid, savename, FILTER=filter, PATH=path, PRINT=print, YRANGE=yrange

  common colori


  IF ~ keyword_set(path) THEN path='./'

;set filename

  obsed=strcompress(path+string(objid)+".obs_sed",/remove_all)
  chidistr=strcompress(path+string(objid)+".pz",/remove_all)
  sed=strcompress(path+string(objid)+".temp_sed",/remove_all)

;read data

  readcol, obsed, lamb_phot, flux_cat, err_flux, temp_flux, skipline=1
  readcol, chidistr, z_val, chi_q, tmpl, skipline=1
  readcol, sed, sed_lambda, sed_flux, skipline=2

  IF keyword_set(FILTER) THEN BEGIN
     


     
  ENDIF 
   


  IF keyword_set(print) THEN m_psopen, savename+"_sed.eps", /maxs, /encapsulate $ 
  ELSE window, 0


  if ~keyword_set(yrange) THEN yrange=[0.,MAX(flux_cat+err_flux)] 
  plot,  lamb_phot, flux_cat, /nodata,yrange=yrange,$
         xtitle=Textoidl("\lambda (A)"), ytitle=Textoidl("F_\nu (normalized)")
  oploterror, lamb_phot, flux_cat, err_flux, psym=5, symsize=3, color=blu, errcolor=blu
  oplot, lamb_phot, temp_flux, psym=6, symsize=3, color=rosso
  oplot, sed_lambda, sed_flux
  
  ;olpot filter if available
  IF keyword_set(FILTER) THEN BEGIN
  ;get file
     readcol, getenv("HOME2")+"/CODES/eazy-1.00/inputs/master.FILTER.RES.info",$
            line, start, finish   
     FOR i=0, N_ELEMENTS(filter)-1 DO BEGIN
        ll=line[where(line EQ filter[i])]
        ss=start[where(line EQ filter[i])] 
        ff=finish[where(line EQ filter[i])]
        ;get actual filter
        readcol, getenv("HOME2")+"/CODES/eazy-1.00/inputs/master.FILTER.RES",$
                 ind_line, lam_filter, trans_filter, skipline=round(ss[0]+1),$
                 numline=round(ss[0]+ff[0])
        ;normalise to maximum
         ;trans_filter=trans_filter*(MAX(flux_cat)/MAX(trans_filter))
        trans_filter=trans_filter*(5D-30/MAX(trans_filter))
        oplot, lam_filter, trans_filter        
     ENDFOR     
  ENDIF 
 

  IF keyword_set(print) THEN m_psclose
  
  IF keyword_set(print) THEN m_psopen, savename+"_chi.eps", /maxs, /encapsulate $ 
  ELSE window, 1
  !y.style=1
  plot, z_val, chi_q, yrange=[MIN(chi_q),MAX(chi_q)],$
        xtitle="z", ytitle=Textoidl("\chi^2"), /ylog
  
  IF keyword_set(print) THEN m_psclose








END

