;+
;
;
;  Calibrate sloan g using gr color
;
;
;
;-


pro calibrate_u_ug, list, airmass, path=path, extinc=extinc, residual=residual, calib=calib, err_calib=err_calib, $
                    save=save, limmag=limmag, outresidual=outresidual


;;calibrate U with U-G color

for i=0, n_elements(list)-1 do begin


    fits=mrdfits(path+list[i],1,/sil)

    ;;transform from ashin to AB mag
    fits.UMAG=asinh2pogs(fits.UMAG,"u")-0.04   ;;correct the U band to true AB

    ;;isolate matched sources (cl = 6 star; cl = 3 galaxy)
    good=where(fits.MATCHM_G ne 0 and fits.MATCHM_U ne 0 and fits.MATCHM_U-fits.MATCHM_G lt 50 and $
               fits.UMAG lt limmag, ngood)

    ;;load good data
    if(i eq 0) then begin
        
        sloan_u=fits[good].UMAG
        sloan_err_u=fits[good].E_UMAG  ; this is imperfect (still ashin) but should be ok or bright star
        
        instr_u=fits[good].MATCHM_U
        instr_err_u=fits[good].ERR_MATCHM_U

        instr_ug=fits[good].MATCHM_U-fits[good].MATCHM_G
        instr_err_ug=sqrt(fits[good].ERR_MATCHM_U^2+fits[good].ERR_MATCHM_G^2.)
        
        airmass_term=replicate(extinc*airmass[i],ngood)
        residual_term=replicate(residual[i],ngood)
        
        class=fits[good].cl

    endif else begin

        sloan_u=[sloan_u,fits[good].UMAG]
        sloan_err_u=[sloan_err_u,fits[good].E_UMAG]

        instr_u=[instr_u,fits[good].MATCHM_U]
        instr_err_u=[instr_err_u,fits[good].ERR_MATCHM_U]

        instr_ug=[instr_ug,fits[good].MATCHM_U-fits[good].MATCHM_G]
        instr_err_ug=[instr_err_ug,sqrt(fits[good].ERR_MATCHM_U^2+fits[good].ERR_MATCHM_G^2.)]

        airmass_term=[airmass_term,replicate(extinc*airmass[i],ngood)]

        residual_term=[residual_term,replicate(residual[i],ngood)]

        class=[class,fits[good].cl]

    endelse
    
endfor


;;do the calibration 
xarray=instr_ug
err_xarray=instr_err_ug
yarray=sloan_u-instr_u-airmass_term-residual_term
err_yarray=sqrt(sloan_err_u^2+instr_err_u^2)

;;plot some 
m_psopen, save, /land
ncl=0


;;restore
new_xarray=xarray
new_yarray=yarray
new_err_yarray=err_yarray
new_err_xarray=err_xarray

;;iter 
for it=0, 5 do begin
    

    ploterror, new_xarray, new_yarray, new_err_xarray, new_err_yarray, psym=1, xtitle='U-V', $
            ytitle=Textoidl('u-U_{ins}-AM'), title='Iter '+rstring(it)

    if(ncl gt 0) then begin
        oplot, [new_xarray[remove]], [new_yarray[remove]], color=fsc_color("green"), psym=1
        new_xarray=new_xarray[clip]
        new_yarray=new_yarray[clip]
        new_err_yarray=new_err_yarray[clip]
        new_err_xarray=new_err_xarray[clip]

    endif
    
    ;calib=LINFIT(new_xarray,new_yarray,CHISQ=chi,prob=prob,MEASURE_ERRORS=new_err_yarray,$
    ;             SIGMA=err_calib,YFIT=yfit) 
    
    calib=ROBUST_LINEFIT(new_xarray,new_yarray,yfit,sig_res,err_calib)
    
    print, "Zero Point ", calib[0], err_calib[0]
    print, "Color Term  ", calib[1], err_calib[1]
    ;print, "Chi  ", chi/(n_elements(new_yarray)-2)
    ;print, "Prob  ", prob

    
    ;;plot fit
    oplot, mkarr(-5,5,0.01), mkarr(-5,5,0.01)*calib[1]+calib[0], line=1, color=fsc_color("red")
    
    xyouts, min(new_xarray), min(new_yarray), $
            "ZP "+string(calib[0],format="(F6.3)")+" +/- "+string(err_calib[0],format="(F6.3)")+$
            "  C "+string(calib[1],format="(F6.3)")+" +/- "+string(err_calib[1],format="(F6.3)")
    

    ;;do clipping for next iter
    if it eq 0 then nsig=10. else  nsig=5.
    clip=where(abs((yfit-new_yarray)/new_err_yarray) lt nsig,ncl,complement=remove) 

endfor

;;replot on original data

    ploterror, xarray, yarray, err_xarray, err_yarray, psym=1, xtitle='U-V', $
         ytitle=Textoidl('u-U_{ins}-AM'), title='Final '+"ZP "+string(calib[0],format="(F6.3)")+$
            " +/- "+string(err_calib[0],format="(F6.3)")+$
            "  C "+string(calib[1],format="(F6.3)")+" +/- "+string(err_calib[1],format="(F6.3)")
    oplot, mkarr(-5,5,0.01), mkarr(-5,5,0.01)*calib[1]+calib[0], line=1, color=fsc_color("red")

    

;;compute the residual for each field

outresidual=residual-residual

for i=0, n_elements(list)-1 do begin
    
    fits=mrdfits(path+list[i],1,/sil)
    fits.UMAG=asinh2pogs(fits.UMAG,"u")-0.04   ;;correct the U band to true AB
    good=where(fits.MATCHM_U ne 0 and fits.MATCHM_G ne 0 and $
               fits.MATCHM_U-fits.MATCHM_G lt 50 and fits.UMAG lt limmag, ngood)

    sloan_u=fits[good].UMAG
    sloan_err_u=fits[good].E_UMAG
    
    instr_u=fits[good].MATCHM_U
    instr_err_u=fits[good].ERR_MATCHM_U
    
    instr_ug=fits[good].MATCHM_U-fits[good].MATCHM_G
    instr_err_ug=sqrt(fits[good].ERR_MATCHM_U^2+fits[good].ERR_MATCHM_G^2.)
    
    airmass_term=replicate(extinc*airmass[i],ngood)

    ucal=instr_u+calib[0]+calib[1]*(instr_ug)+airmass_term+residual[i]
    ucal_err=sqrt(instr_err_u^2+err_calib[0]^2+instr_err_ug^2)
    res=sloan_u-ucal
    res_err=(ucal_err^2+sloan_err_u^2)
    
    outresidual[i]=median(res)
    
    
    ploterror,  sloan_u, res, sloan_err_u, res_err,  psym=1, title=list[i]+" Res: "+string(median(res),format='(F5.2)'), $
            xtitle='sloan u', ytitle='Residual', yrange=[median(res)-5*stddev(res),median(res)+5*stddev(res)]
    
    oplot, [-100,100], [median(res),median(res)], line=3
    oplot, [-100,100], [median(res),median(res)]+stddev(res), line=1
    oplot, [-100,100], [median(res),median(res)]-stddev(res), line=1

    print, "Residual ", list[i], median(res), stddev(res)

endfor


;;plot the distribution of residuals
if(n_elements(list) gt 1) then begin 
    plothist, outresidual, bin=0.05, xtitle='Residuals'
    print, "Deviations of residual ", stddev(outresidual)
endif

m_psclose



end
