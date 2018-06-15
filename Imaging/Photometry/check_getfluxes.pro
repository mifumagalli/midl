;+
;
;
;  This is a procedure that makes some plot to check the photometry
;  extraction
;
;   objname   -->  name of an object structure
;   path      -->  location of the object stracture
;   sdss      -->  set this if there are sdss info to make 
;                  additional plot
;
;-

pro check_getfluxes, objname, sdss=sdss, path=path

;;get data
if not keyword_set(path) then path='./'
obj=mrdfits(path+objname,1,/sil)

;;open plot
plname='phot_'+(strmid(objname,0,strpos(objname,'.fits')))+'.ps'
m_psopen, plname, /land

;;loop over filters
nfilt=n_elements(obj[0].number)
nobj=n_elements(obj.number[0])



for f=0, nfilt-1 do begin

    ;;check iso mag
    good=where(obj.COLOR_SN[f] gt 0,ngoo)
    if(ngoo gt 0.) then $
      plothist, mk_finite(obj[good].COLOR_FLUX[f]/obj[good].FLUX_ISO[f]), bin=0.05,$
      xrange=[-2,3], title='Color flux for filter '+obj[0].FILTER[f], $
      xtitle='COLOR/ISO', ytitle='Over '+rstring(ngoo)
    
    ;;check total flux
    good=where(obj.TOT_SN[f] gt 0,ngoo)
    if(ngoo gt 0.) then $
      plothist, mk_finite(obj[good].TOT_FLUX_COR[f]/obj[good].FLUX_AUTO[f]), bin=0.05,$
      xrange=[0,3], title='Total flux for filter '+obj[0].FILTER[f], $
      xtitle='TOT_CORR/AUTO', ytitle='Over '+rstring(ngoo)
    
    ;;check SN and upper limits
    plot, obj.TOT_FLUX_COR[f], obj.TOT_SN[f], /xlog, xrange=[1,max(obj.TOT_FLUX_COR[f])], $
      /ylog, yrange=[0.1,max(obj.TOT_SN[f])], psym=3, $
      title='S/N for filter '+obj[0].FILTER[f], xtitle='TOT_CORR', $
      ytitle='S/N'
    
    if(f gt 0) then begin
        ;;make color color comparison
        good=where(obj.COLOR_SN[f] gt 0 and obj.COLOR_SN[0] gt 0. and $
                   obj.TOT_SN[f] gt 0 and obj.TOT_SN[0] gt 0.,ngoo)
        
        if(ngoo gt 0) then begin
            plot, obj[good].TOT_FLUX_COR[f]/obj[good].TOT_FLUX_COR[0], $
              obj[good].COLOR_FLUX[f]/obj[good].COLOR_FLUX[0], $
              title='Color '+obj[0].FILTER[f]+' - '+obj[0].FILTER[0], psym=1, $
              xtitle='TOT_CORR', ytitle='ISO', xrange=[0,10], yrange=[0,10]
            oplot, [-100,100], [-100,100], line=2
        endif
    endif
    
    ;;make a SDSS comparison
    if keyword_set(sdss) then begin
        
        ;;position
        good=where(obj.SDSS_RAJ2000[f] gt 0.)
        
        plot, obj[good].SDSS_RAJ2000[f], obj[good].SDSS_DEJ2000[f], $
          psym=1, xtitle='RA', ytitle='DEC', /ynozero, $
          title='SDSS position for filter '+obj[0].FILTER[f]
        oplot, obj[good].ALPHA_J2000[f], obj[good].DELTA_J2000[f], $
          psym=4, color=250
        
        ;;find residual
        res=djs_diff_angle(obj[good].SDSS_RAJ2000[f],obj[good].SDSS_DEJ2000[f],$
                           obj[good].ALPHA_J2000[f], obj[good].DELTA_J2000[f])
        
        medres=median(res*3600.)
        
        plothist, res*3600., bin=0.01, $
          title='Match position for filter '+obj[0].FILTER[f], $
          xtitle='Residual in a.s.'
        
        oplot, [medres,medres], [-1,10000], line=2

    endif
    
endfor

m_psclose

end

