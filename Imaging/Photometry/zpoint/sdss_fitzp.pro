;+
;  Photometric calibrate an object structure by comparing against sdss
;  Result is an AB mag.
;   
;  If the keyword preserve is set, the code uses stars with neutral 
;  color to set a ZP for the filter in used, otherwise the final magnitudes
;  are calibrated over the ugriz system.  
;  Also, a correction for galactic extinction is applied 
;
;  objname   --> name of a file structure 
;  path      --> where to find the object structure
;  instr     --> the instrument used (LRIS,LBC supported) to load  
;               defualt values of AM.
;  am        --> an array of air masses coeficient (1 per filter) to be used in
;               the fit. If not provided, default are used. 
;  obs_am    --> an array (1 per filter) of air masses at the moment of the observation
;  sdss_fil  --> an array of sdss filer used during the sdss query                            
;  sdss_maglim --> if set to a 2D array, it allows you to tweak the
;                  limit in magnitude used for the ZP fit. Each
;                  element of the array is the faint and bright limit.
;                  E.g. for different filters [[22,17],[21,15],[23,14]]
;
;  Currently, it is implemented only a calibration with fixed AM
;  and free color term. This is far from being general.. One could add calibrations
;  without AM or with free AM and color terms.
;
;  To be improved:
;     -- Check the optimal range for ZP fitting (look at http://adsabs.harvard.edu/abs/2008AJ....135..264C)
;   
;
;-


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;This is the utility function that shifts over sdss system
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

pro sdss_bootstrap, obj, f, nfilt, sdss_fil, obs_am, am, lun, nobj, sdss_maglim=sdss_maglim


printf, lun, 'Shifting over sdss system'

;;define index color
if(f lt nfilt-1) then begin
    c2=f+1
    c1=f
endif else begin
    c2=f
    c1=f-1
endelse

splog, 'Calibrate ', obj[0].filter[f], ' with color ', obj[0].filter[c1], '-',  obj[0].filter[c2]
printf, lun, '-------------------------------------'
printf, lun, 'Calibrate ', obj[0].filter[f], ' with color ', obj[0].filter[c1], '-',  obj[0].filter[c2]

;;isolate object that are in sdss and with positive total flux in
;;the filter and in the colors. Limit to sdss mag < 22. Exclude
;;too bright, possibly saturated and with some PSF modeling problems
if not keyword_set(sdss_maglim) then begin
    faint_mag=22
    bright_mag=16
endif else begin
    faint_mag=sdss_maglim[0,f]
    bright_mag=sdss_maglim[1,f]
endelse



obj_cal=where(obj.sdss_mag[f] gt bright_mag and obj.sdss_mag[f] le faint_mag and $
              obj.tot_flux_cor[f] gt 0. and obj.color_flux[c1] gt 0. and obj.color_flux[c2] gt 0.,nus)

splog, 'Found ', nus, ' objects for calibrations'

;;extract magnitude
sdss_mag_ash=obj[obj_cal].sdss_mag[f]
obj_mag=-2.5*alog10(obj[obj_cal].tot_flux_cor[f])
obj_c1=-2.5*alog10(obj[obj_cal].color_flux[c1])
obj_c2=-2.5*alog10(obj[obj_cal].color_flux[c2])
obj_color=(obj_c1-obj_c2)
;;this is an approximation of a asym err bar
obj_errmag=1.0857/obj[obj_cal].tot_sn[f]


;;find linear magnitude and clip negative fluxes
sdss_mag=asinh2pogs(sdss_mag_ash,sdss_fil[f],/clip,used=used)
;;apply correction form SDSS magnitude to AB magnitude
if(sdss_fil[f] eq 'u') then sdss_mag=sdss_mag-0.04
if(sdss_fil[f] eq 'z') then sdss_mag=sdss_mag+0.02

obj_mag=obj_mag[used]
obj_c1=obj_c1[used]
obj_c2=obj_c2[used]
obj_color=obj_color[used]
obj_errmag=obj_errmag[used]


;;fit zp with color term and fixed am
;;This can be iterated to improve the solution
;;nb x_photsol2 applies Am term in mag_instr
x_photsol2, obj_mag, sdss_mag, obj_errmag, obs_am[f], $
  coeffs, sigma_coeff, setam=am[f], color=obj_color,  CHISQ=chi 


splog, 'ZP: ', coeffs[0]
splog, 'eZP: ', sigma_coeff[0]
splog, 'AM: ', am[f]
splog, 'eAM: ', 0.
splog, 'CLR: ', coeffs[1]
splog, 'eCLR: ', sigma_coeff[1]
splog, 'Chi^2/dof: ', chi/(nus-2) 

printf, lun, 'ZP: ', coeffs[0]
printf, lun, 'eZP: ', sigma_coeff[0]
printf, lun, 'AM: ', am[f]
printf, lun, 'eAM: ', 0.
printf, lun, 'CLR: ', coeffs[1]
printf, lun, 'eCLR: ', sigma_coeff[1]
printf, lun, 'Chi^2/dof: ', chi/(nus-2) 

;;nb x_photsol2 already includes AM term in mag_instr
obj_magcal=obj_mag+coeffs[0]-coeffs[1]*obj_color
residual=obj_magcal-sdss_mag
djs_iterstat, residual, mean=all_res, mask=mask
splog, 'Residual ', all_res
printf, lun, 'Residual ', all_res


;;plot calibrated
plot, sdss_mag, obj_magcal, psym=1, /ynozero, title='Calibration for '+obj[0].filter[f], $
  xtitle='SDSS MAG', ytitle='CALIBRATED MAG'
diag=mkarr(0,30,0.1)
oplot, diag, diag, line=1

;;plot residual
bright_r=residual[where(sdss_mag lt 20)]
bright_m=sdss_mag[where(sdss_mag lt 20)]
djs_iterstat, bright_r, mean=bri_res
splog, 'Residual bright ', bri_res
printf, lun, 'Residual bright ', bri_res

plot, sdss_mag, residual, psym=1, /ynozero, title='Residual '+rstring(all_res), $
  xtitle='sdss mag', ytitle='all residual'
plot, bright_m, bright_r, psym=1, /ynozero, yrange=[-1,1], title='Residual '+rstring(bri_res), $
  xtitle='sdss mag', ytitle='bright residual'

;;plot color term
plot, sdss_mag, coeffs[1]*obj_color, psym=1, /ynozero, $
  title=string(obj[0].filter[c1],'-',obj[0].filter[c2]), $
  xtitle='sdss mag', ytitle='Color term'

;;find galactic extinction based on position at the new frequency
;;used (sdss filters)
;;fix filter name
if(sdss_fil[f] eq 'u') then  ext_filter='usdss'
if(sdss_fil[f] eq 'g') then  ext_filter='gsdss'
if(sdss_fil[f] eq 'r') then  ext_filter='rsdss'
if(sdss_fil[f] eq 'i') then  ext_filter='isdss'
if(sdss_fil[f] eq 'z') then  ext_filter='zsdss'

ebv=ebv_dust(obj.alpha_j2000[0],obj.delta_j2000[0],a=a_ext,band=ext_filter,/gg)

;;store all the calibration
obj.galact_a[f]=a_ext
obj.zp_mag[f]=replicate(coeffs[0],nobj)
;;If color is lt 0, suppress color term in calibration
obj.color_term[f]=mk_finite(coeffs[1]*(-2.5*alog10(obj.color_flux[c1])+2.5*alog10(obj.color_flux[c2])))
obj.am_term[f]=am[f]*obs_am[f]

;;compute final mag, calibrated and apply reddening correction
;;If the object has a SN=0, i.e. zero or negative flux, put 
;;the magnitude to 99 and set the error bar to 1sigma limiting magntiude
;;If the object is not measured (marked as defect by sextractor, set
;;error to zero)

;;First object with flux
detect_color=where(obj.color_sn[f] gt 0 and obj.color_flux[f] gt 0)
detect_total=where(obj.tot_sn[f] gt 0 and obj.tot_flux_cor[f] gt 0)

obj[detect_total].tot_magcal_cor[f]=-2.5*alog10(obj[detect_total].tot_flux_cor[f])+$
  obj[detect_total].zp_mag[f]-obj[detect_total].color_term[f]-obj[detect_total].am_term[f]$
  -obj[detect_total].galact_a[f]

obj[detect_color].color_magcal_cor[f]=-2.5*alog10(obj[detect_color].color_flux[f])+$
  obj[detect_color].zp_mag[f]-obj[detect_color].color_term[f]-obj[detect_color].am_term[f]$
  -obj[detect_color].galact_a[f]

;;compute error adding rms calibration 
obj[detect_total].errtot_magcal_cor[f]=sqrt((1.0857/obj[detect_total].tot_sn[f])^2+all_res^2)
obj[detect_color].errcolor_magcal_cor[f]=sqrt((1.0857/obj[detect_color].color_sn[f])^2+all_res^2)

;;now set the non detection and upper limits for zero negative flux 
;;Skip non measured with zero flux
nodet_color=where(obj.color_sn[f] le 0 and obj.color_flux[f] lt 0)
nodet_total=where(obj.tot_sn[f] le 0 and obj.tot_flux_cor[f] lt 0)

obj[nodet_total].tot_magcal_cor[f]=99
obj[nodet_color].color_magcal_cor[f]=99

;;set the error bar
obj[nodet_total].errtot_magcal_cor[f]=sqrt((-2.5*alog10(obj[nodet_total].errtot_flux_cor[f])+$
                                                obj[nodet_total].zp_mag[f]-obj[nodet_total].color_term[f]$
                                                -obj[nodet_total].am_term[f]$
                                                -obj[nodet_total].galact_a[f])^2+all_res^2)

obj[nodet_color].errcolor_magcal_cor[f]=sqrt((-2.5*alog10(obj[nodet_color].errcolor_flux[f])+$
                                                  obj[nodet_color].zp_mag[f]-obj[nodet_color].color_term[f]$
                                                  -obj[nodet_color].am_term[f]$
                                                  -obj[nodet_color].galact_a[f])^2+all_res^2)

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;Possible combination:
;;
;;SN > 0 both mag and err have a value
;;SN < 0   --> never observed since marked as defect by sext has 0 in
;;             both mag and error 
;;         --> undetected has 99 in magnitude and upper limit in error
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

end

pro sdss_preserve, obj, f, nfilt, sdss_fil, obs_am, am, lun, nobj, instr, sdss_maglim=sdss_maglim

printf, lun, 'Keep original photmetric system'

;;define index color in sdss
if(f lt nfilt-1) then begin
    c2=f+1
    c1=f
endif else begin
    c2=f
    c1=f-1
endelse
;;fix if same color 
if(sdss_fil[c1] eq sdss_fil[c2] and c2 lt nfilt-1) then c2++  
if(sdss_fil[c1] eq sdss_fil[c2] and c2 ge nfilt-1) then c2--  

splog, 'Calibrate ', obj[0].filter[f], ' for 0 color ', sdss_fil[c1], '-', sdss_fil[c2]
printf, lun, '-------------------------------------'
printf, lun, 'Calibrate ', obj[0].filter[f], ' for 0 color ', sdss_fil[c1], '-', sdss_fil[c2]

;;isolate object that are in sdss and with positive total flux in
;;the filter and in the colors. Limit to sdss mag < 22. Exclude
;;too bright, possibly saturated and with some PSF modeling problems
;;Be a bit more specific with u band
if not keyword_set(sdss_maglim) then begin
    if(sdss_fil[f] eq 'u') then begin
        faint_mag=22
        bright_mag=14
    endif else begin
        faint_mag=22
        bright_mag=16
    endelse
endif else begin
    faint_mag=sdss_maglim[0,f]
    bright_mag=sdss_maglim[1,f]
    faint_mag_c1=sdss_maglim[0,c1]
    bright_mag_c1=sdss_maglim[1,c1]
    faint_mag_c2=sdss_maglim[0,c2]
    bright_mag_c2=sdss_maglim[1,c2]
 
    splog, 'SDSS mag. limit ', faint_mag, bright_mag
endelse

;obj_cal=where(obj.sdss_mag[f] gt bright_mag and obj.sdss_mag[f] le faint_mag and $
;              obj.tot_flux_cor[f] gt 0. and obj.sdss_mag[c1] gt bright_mag_c1 and obj.sdss_mag[c1] le faint_mag_c1 and $
;              obj.sdss_mag[c2] gt bright_mag_c2 and obj.sdss_mag[c2] le faint_mag_c2, nus)

obj_cal=where(obj.sdss_mag[f] gt bright_mag and obj.sdss_mag[f] le faint_mag and $
              obj.tot_flux_cor[f] gt 0., nus)


splog, 'Found ', nus, ' objects for calibrations'

;;extract magnitude
sdss_mag_ash=obj[obj_cal].sdss_mag[f]
obj_mag=-2.5*alog10(obj[obj_cal].tot_flux_cor[f])
sdss_c1_ash=obj[obj_cal].sdss_mag[c1]
sdss_c2_ash=obj[obj_cal].sdss_mag[c2]
;;this is an approximation of a asym err bar
obj_errmag=1.0857/obj[obj_cal].tot_sn[f]

;;find linear magnitude and clip negativre fluxes
sdss_mag=mk_finite(asinh2pogs(sdss_mag_ash,sdss_fil[f]))
sdss_c1=mk_finite(asinh2pogs(sdss_c1_ash,sdss_fil[c1]))
sdss_c2=mk_finite(asinh2pogs(sdss_c2_ash,sdss_fil[c2]))

used=where(sdss_mag gt 0 and sdss_c1 gt 0 and sdss_c2 gt 0)

;;clip
obj_mag=obj_mag[used]
obj_errmag=obj_errmag[used]
sdss_mag=sdss_mag[used]
sdss_c1=sdss_c1[used]
sdss_c2=sdss_c2[used]


;;apply correction form SDSS magnitude to AB magnitude
if(sdss_fil[f] eq 'u') then sdss_mag=sdss_mag-0.04
if(sdss_fil[f] eq 'z') then sdss_mag=sdss_mag+0.02
if(sdss_fil[c1] eq 'u') then sdss_c1=sdss_c1-0.04
if(sdss_fil[c1] eq 'z') then sdss_c1=sdss_c1+0.02
if(sdss_fil[c2] eq 'u') then sdss_c2=sdss_c2-0.04
if(sdss_fil[c2] eq 'z') then sdss_c2=sdss_c2+0.02

;;define color
sdss_color=(sdss_c1-sdss_c2)

;;find all zp compensating for AM 
all_zp=sdss_mag-obj_mag+am[f]*obs_am[f]

;;Fit a zp / color linear relation with robust fit.
;;Your ZP is where the line goes through zero
coeff=robust_linefit(sdss_color,all_zp,model,std_res,sigma_coeff)
;;clean outliers and do second fit
qdone=djs_reject(all_zp,model,outmask=out,upper=4.,lower=4.)
good=where(out gt 0)
coeff=robust_linefit(sdss_color[good],all_zp[good],model,std_res,sigma_coeff)


splog, 'ZP: ', coeff[0]
splog, 'eZP: ', sigma_coeff[0]
splog, 'AM: ', am[f]
splog, 'eAM: ', 0.

printf, lun, 'ZP: ', coeff[0]
printf, lun, 'eZP: ', sigma_coeff[0]
printf, lun, 'AM: ', am[f]
printf, lun, 'eAM: ', 0.

;;make some plot of the fit
plot, sdss_color[good], all_zp[good], psym=1, /ynozero, $
  title='Calibration for '+obj[0].filter[f]+' ZP '+rstring(coeff[0]), $
  xtitle='Color '+sdss_fil[c1]+'-'+sdss_fil[c2], ytitle='ZP-AM'
oplot, sdss_color[good], model


;;find galactic extinction for the riginal filter
;;fix lbc names
ext_filter=obj[0].filter[f]
if(instr eq 'LBC') then begin
    if(ext_filter eq 'SD') then  ext_filter='U'
    if(ext_filter eq 'B-') then  ext_filter='B'
    if(ext_filter eq 'V-') then  ext_filter='V'
    if(ext_filter eq 'R-') then  ext_filter='R'
    if(ext_filter eq 'I-') then  ext_filter='I'
endif

ebv=ebv_dust(obj.alpha_j2000[0],obj.delta_j2000[0],a=a_ext,band=ext_filter,/gg)

;;store all the calibration
obj.galact_a[f]=a_ext
obj.zp_mag[f]=replicate(coeff[0],nobj)
;;color term is zero
obj.color_term[f]=obj.galact_a[f]-obj.galact_a[f]
obj.am_term[f]=am[f]*obs_am[f]


;;compute final mag, calibrated and apply reddening correction
;;If the object has a SN=0, i.e. zero or negative flux, put 
;;the magnitude to 99 and set the error bar to 1sigma limiting magntiude
;;If the object is not measured (marked as defect by sextractor, set
;;error to zero)

;;First object with positive flux
detect_color=where(obj.color_sn[f] gt 0 and obj.color_flux[f] gt 0)
detect_total=where(obj.tot_sn[f] gt 0 and obj.tot_flux_cor[f] gt 0)

obj[detect_total].tot_magcal_cor[f]=-2.5*alog10(obj[detect_total].tot_flux_cor[f])+$
  obj[detect_total].zp_mag[f]-obj[detect_total].color_term[f]-obj[detect_total].am_term[f]$
  -obj[detect_total].galact_a[f]

obj[detect_color].color_magcal_cor[f]=-2.5*alog10(obj[detect_color].color_flux[f])+$
  obj[detect_color].zp_mag[f]-obj[detect_color].color_term[f]-obj[detect_color].am_term[f]$
  -obj[detect_color].galact_a[f]

;;compute error adding rms calibration 
obj[detect_total].errtot_magcal_cor[f]=sqrt((1.0857/obj[detect_total].tot_sn[f])^2+sigma_coeff[0]^2)
obj[detect_color].errcolor_magcal_cor[f]=sqrt((1.0857/obj[detect_color].color_sn[f])^2+sigma_coeff[0]^2)

;;now set the non detection and upper limits for zero negative flux 
;;Skip those with identical zero flux which are not measured
nodet_color=where(obj.color_sn[f] le 0 and obj.color_flux[f] lt 0)
nodet_total=where(obj.tot_sn[f] le 0 and obj.tot_flux_cor[f] lt 0)

obj[nodet_total].tot_magcal_cor[f]=99
obj[nodet_color].color_magcal_cor[f]=99

;;set error bar
obj[nodet_total].errtot_magcal_cor[f]=sqrt((-2.5*alog10(obj[nodet_total].errtot_flux_cor[f])+$
                                                obj[nodet_total].zp_mag[f]-obj[nodet_total].color_term[f]$
                                                -obj[nodet_total].am_term[f]$
                                                -obj[nodet_total].galact_a[f])^2+sigma_coeff[0]^2)

obj[nodet_color].errcolor_magcal_cor[f]=sqrt((-2.5*alog10(obj[nodet_color].errcolor_flux[f])+$
                                                  obj[nodet_color].zp_mag[f]-obj[nodet_color].color_term[f]$
                                                  -obj[nodet_color].am_term[f]$
                                                  -obj[nodet_color].galact_a[f])^2+sigma_coeff[0]^2)


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;Possible combination:
;;
;;SN > 0 both mag and err have a value
;;SN < 0   --> never observed since marked as defect by sext has 0 in
;;             both mag and error 
;;         --> undetected has 99 in magnitude and upper limit in error
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;



;;plot the magnitude distribution
plothist, obj.tot_magcal_cor[f], bin=0.1, xrange=[15,30], xtitle='Cal Mag', $
  title='Calibration for '+obj[0].filter[f]

;;plot a comparison with sdss
ind=where(obj.tot_sn[f] gt 0 and obj.sdss_mag[f] gt 0, ncfr)
if(ncfr gt 0) then begin
    plot, obj[ind].sdss_mag[f], obj[ind].tot_magcal_cor[f]-obj[ind].sdss_mag[f], psym=1, $
      title='Calibration for '+obj[0].filter[f], xtitle=sdss_fil[f], ytitle=string(obj[0].filter[f],'-',sdss_fil[f])
    oploterror, obj[ind].sdss_mag[f], obj[ind].tot_magcal_cor[f]-obj[ind].sdss_mag[f], obj[ind].errtot_magcal_cor[f], psym=3
    oplot, [-100,100], [0,0], line=2
    
endif
end

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;This is the main procedure
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

pro sdss_fitzp, objname, path=path, am=am, obs_am=obs_am, instr=instr, $
                sdss_fil=sdss_fil, preserve=preserve, sdss_maglim=sdss_maglim


;;default
if not keyword_set(path) then path='./'
if not keyword_set(instr) then instr='LRIS'
if not keyword_set(sdss_fil) then begin
    sdss_fil=['u','g','r','i']
    splog, 'Assuming the following filters...'
endif

;;open data structure
obj=mrdfits(path+objname,1,/sil)
nfilt=n_elements(obj[0].number)
nobj=n_elements(obj.number[0])
obs_filt=obj[0].filter


;;clean any previous calibration
obj.tot_magcal_cor=0.
obj.color_magcal_cor=0.
obj.errtot_magcal_cor=0.
obj.errcolor_magcal_cor=0.
obj.galact_a=0.
obj.zp_mag=0.
obj.color_term=0.
obj.am_term=0.


;;load default am
if not keyword_set(obs_am) then obs_am=replicate(1.,nfilt)
if not keyword_set(am) then begin
    
    am=replicate(1.,nfilt)
    splog, 'Loading default AM'
    
    ;;open inst specific list
    if(instr eq 'LRIS') then $
      readcol, getenv("MIDL")+'/Imaging/Photometry/utility/MaunaKea_am.txt', $
      am_filt, am_val,format='A,F', /sil
    ;;lbt
    if(instr eq 'LBC') then $
      readcol, getenv("MIDL")+'/Imaging/Photometry/utility/LBT_am.txt', $
      am_filt, am_val,format='A,F', /sil
    for kk=0, nfilt-1 do begin
        ;;extract am
        fndfilt=where(am_filt eq rstring(obs_filt[kk]))
        if(fndfilt gt -1) then  am[kk]=am_val[fndfilt]
    endfor
   
    splog, 'Using AM ', am
    
endif

;;open ps
m_psopen, 'zpcal_'+strmid(objname,0,strpos(objname,'.fits'))+'.ps', /land
;;open txt
openw, lun, 'zpcal_'+strmid(objname,0,strpos(objname,'.fits'))+'.txt', /get_lun


;;loop over filter
for f=0, nfilt-1 do begin
    ;;compute the calibration
    if keyword_set(preserve) then sdss_preserve, obj, f, nfilt, sdss_fil, $
      obs_am, am, lun, nobj, instr, sdss_maglim=sdss_maglim $
    else sdss_bootstrap, obj, f, nfilt, sdss_fil, obs_am, am, lun, nobj, $
      sdss_maglim=sdss_maglim
endfor 

;;close ps and txt
m_psclose
close, /all

;;save structure
mwrfits, obj, path+objname, /cre

end
