;+
;
;  Given a pair of reduced pixel tables (object and sky exposures), produces
;  a SKYCORR scaled sky model for each IFU.  
;
;  pixtabo   -> name of the object pixel table (with illum corr)
;  pixtabs   -> name of the sky pixel table (with illum corr)
;  psky     -> percentiles of pixels used to extract the sky spectrum in the 
;              object frame  [min/max]. In the sky frame we always use [45-55].
;  quick    -> select a narrow wave range to perform parameter study 
;
;-


function make_skyspec, data, var, lambda, nstep, npix, lambin, skysampl, wminmax, psky
    
    modelsky_wave=mkarr(wminmax[0]-1,wminmax[1]+1,skysampl)
    modelsky_flux=fltarr(n_elements(modelsky_wave),nstep)
    
    ;;find blocks of wave - do it by index which is way faster
    nbloc=floor((wminmax[1]-wminmax[0])/lambin)
    indxbloc=floor(npix/nbloc)
    offindx=1.*indxbloc/nstep

    ;;run the loop - needed to avoid 'edges' effect do to skylines 
    for steps=0, nstep-1 do begin
       
       splog, 'Iteration ', steps+1, ' out of ',  nstep
       
       ;;create storage for pixels
       skypixflux=0.
       skypixwave=0.
       ;skypixvar=0.

       ;;loop over blocks 
       for bb=0, nbloc-2 do begin

    	  b_data=data[bb*indxbloc+offindx*steps:(bb+1)*indxbloc-1+offindx*steps]
    	  b_var=var[bb*indxbloc+offindx*steps:(bb+1)*indxbloc-1+offindx*steps]
    	  b_lambda=lambda[bb*indxbloc+offindx*steps:(bb+1)*indxbloc-1+offindx*steps]
    	  
    	  ;;find percentiles
    	  p=percentiles(b_data,value=psky)
    	  inside=where(b_data ge p[0] and b_data lt p[1],nskyp)
    	  skypixflux=[skypixflux,b_data[inside]]
    	  skypixwave=[skypixwave,b_lambda[inside]]
    	  ;skypixvar=[skypixvar,b_var[inside]]
    			
       endfor

       splog, 'Construct sky spectrum for this step...'

       ;;construct a regularized version of the sky
       modelsky_flux_tmp=(modelsky_flux[*,0])*0.
       modelsky_num_tmp=(modelsky_flux[*,0])*0.
       
       skypixindex=(skypixwave-wminmax[0]+1)/skysampl
       populate_image, modelsky_flux_tmp, skypixindex, weights=skypixflux
       populate_image, modelsky_num_tmp, skypixindex
       modelsky_flux[*,steps]=modelsky_flux_tmp/modelsky_num_tmp

       ;;check 
       plot, skypixwave, skypixflux, psym=1, xrange=[6790,6890]
       oplot, modelsky_wave, modelsky_flux[*,steps], line=1, color=fsc_color('red')

    endfor
  
    splog, 'Construct master spectrum...'
  
    return, median(modelsky_flux,dim=2)
 
 end   

pro make_skycorrpfile, dir, datahead, ifuslice=ifuslice

openw, lun, dir+'/skycorr.par', /get_lun
printf, lun, 'INST_DIR='+getenv('SKYCORR_DIR') ;;/usr/local/software/skycorr/current' ;To be replaced or use an env variable?
printf, lun, 'INPUT_OBJECT_SPECTRUM='+dir+'/mastersky_obj.txt'
printf, lun, 'INPUT_SKY_SPECTRUM='+dir+'/mastersky_sky.txt'
printf, lun, 'OUTPUT_DIR='+dir
printf, lun, 'OUTPUT_NAME=skycorr'    ;_ifu'+ifuslice
printf, lun, 'COL_NAMES=lambda flux NONE NONE'
printf, lun, 'DEFAULT_ERROR=0.01'
printf, lun, 'WLG_TO_MICRON=1.e-4'
printf, lun, 'VAC_AIR=air'
printf, lun, 'DATE_VAL='+string(fxpar_hier(datahead, 'MJD-OBS'))
printf, lun, 'TIME_VAL='+string(fxpar_hier(datahead, 'UTC'))
printf, lun, 'TELALT_VAL='+string(fxpar_hier(datahead, 'TEL.ALT', /hier))
printf, lun, 'LINETABNAME=airglow_groups.dat'
printf, lun, 'VARDATNAME=airglow_var.dat'
printf, lun, 'SOLDATURL=ftp.geolab.nrcan.gc.ca/data/solar_flux/monthly_averages'
printf, lun, 'SOLDATNAME=solflux_monthly_average.txt'
printf, lun, 'SOLFLUX=-1'
printf, lun, 'FWHM=5.0'
printf, lun, 'VARFWHM=0'
printf, lun, 'LTOL=1e-2'
printf, lun, 'MIN_LINE_DIST=2.5'
printf, lun, 'FLUXLIM=-1'
printf, lun, 'FTOL=1e-3'
printf, lun, 'XTOL=1e-3'
printf, lun, 'WTOL=1e-3'
printf, lun, 'CHEBY_MAX=7'
printf, lun, 'CHEBY_MIN=3'
printf, lun, 'CHEBY_CONST=0.'
printf, lun, 'REBINTYPE=1'
printf, lun, 'WEIGHTLIM=0.67'
printf, lun, 'SIGLIM=15'
printf, lun, 'FITLIM=0.'
printf, lun, 'PLOT_TYPE=N'
free_lun, lun

end

pro muse_skytweak, pixtabo, pixtabs, psky=psky, quick=quick

  ;;default 
  if ~keyword_set(psky) then psky=[0.20,0.30]
  
  ;;load pixel table 
  muse_loadpixtab, pixtabo, outhead=headobj, outstr=strobj
  muse_loadpixtab, pixtabs, outhead=headsky, outstr=strsky
  
  ;;loop over steps to fill a window
  lambin=0.4D ;;do chunks of 0.4A
  nstep=4     ;;cover the full range in 4 steps of 0.1A
  skysampl=0.4  ;;0.4 is optimal sampling in pixel table - careful with changes 


  if keyword_set(quick) then begin

     splog, 'Trim data...'
     sel=where(strobj.lambda gt 6200 and strobj.lambda lt 6900)
     strobj.lambda=strobj.lambda[sel]
     strobj.data=strobj.data[sel]
     strobj.dq=strobj.dq[sel]
     strobj.stat=strobj.stat[sel]
     
     str={lambda:strobj.lambda[sel],data:strobj.data[sel],dq:strobj.dq[sel],stat:strobj.stat[sel],$
          xpix:strobj.xpix[sel],ypix:strobj.ypix[sel],origin:strobj.origin[sel]}
     
  endif
  
  ;;find min max in object and sky (check if the instrument mode is consistent)
  wminmaxo=minmax(strobj.lambda)
  wminmaxs=minmax(strsky.lambda)
  
  if (wminmaxo[0] ne wminmaxs[0] or wminmaxo[1] ne wminmaxs[1]) then begin
   splog, 'The lambda range of the object and sky cubes does not match. Is the instrument mode the same?'
   stop
  endif
  
  
  ;for ifu =1,24 do begin
    
  ;  splog, 'Now processing IFU...'+string(ifu)
  ;  for slice=1,48 do begin
    
    ;;sort the good data to speed up loop by index 
    splog, 'Creating master object and sky spectra'
    
    ;goodpix=where(strobj.dq eq 0 and strobj.ifu eq ifu and str.slice eq slice,npixo)
    goodpix=where(strobj.dq eq 0 ,npixo)
    data=strobj.data[goodpix]
    var=strobj.stat[goodpix]
    lambda=strobj.lambda[goodpix]
  
    lsort=sort(lambda)
    data=data[lsort]
    var=var[lsort]
    lambda=lambda[lsort]
    
    ;;clean a bit 
    undefine, lsort
    
    ;goodpixs=where(strsky.dq eq 0 and strsky.ifu eq ifu and str.slice eq slice,npixs)
    goodpixs=where(strsky.dq eq 0 ,npixs)
    datas=strsky.data[goodpix]
    vars=strsky.stat[goodpix]
    lambdas=strsky.lambda[goodpix]
  
    lsort=sort(lambdas)
    datas=datas[lsort]
    vars=vars[lsort]
    lambdas=lambdas[lsort]
    
    ;;clean a bit 
    undefine, lsort

    ;;make space for the sky models
    
    ;;run the function - needed to avoid 'edges' effect do to skylines 
    masterobj_flux =  make_skyspec( data,  var,  lambda,  nstep, npixo, lambin, skysampl, wminmaxo, psky)
    mastersky_flux =  make_skyspec( datas, vars, lambdas, nstep, npixs, lambin, skysampl, wminmaxo, [0.45,0.55])
    mastersky_wave =  mkarr(wminmaxo[0]-1,wminmaxo[1]+1,skysampl)
 
    ;Here we need to call SKYCORR
    
    cd, current=cwd
    ifudir = cwd ;+'/IFU'+strtrim(string(ifu,format='(I02)'),2)+'_slice'+strtrim(string(slice, format='(I02)'),2)
    ;ifuslice = strtrim(string(ifu,format='(I02)'),2)+'_'+strtrim(string(slice, format='(I02)'),2)
    ;file_mkdir, ifudir
    
    ;write object and sky spectra
    splog, 'Running skycorr...'
    
    savenameo=ifudir+'/mastersky_obj.txt'
    savenames=ifudir+'/mastersky_sky.txt'
    forprint, mastersky_wave, masterobj_flux, textout=savenameo, comment='#LAMBDA FLUX'
    forprint, mastersky_wave, mastersky_flux, textout=savenames, comment='#LAMBDA FLUX'
    
    ;prepare the SKYCORR par file and run it.
    make_skycorrpfile, ifudir, headsky.headmain
    spawn, 'skycorr skycorr.par'
    
    ;Read the results of the subtraction
    readcol, ifudir+'/mastersky_obj_SC.txt', lamcor, fluxcor, fluxSCcor, comment='#'
    
    corrsky_flux = fluxcor-fluxSCcor
    
    ;check 
    plot, mastersky_wave, mastersky_flux, xrange=[6790,6890]
    oplot,mastersky_wave, masterobj_flux, color=fsc_color('red')
  
    splog, 'Interpolate sky over pixel table...'
    skysubval=interpol(corrsky_flux,mastersky_wave,strobj.lambda[goodpix],/spline)
  
    ;subtract
    strobj.data[goodpix]=strobj.data[goodpix]-skysubval
    strobj.stat[goodpix]=strobj.stat[goodpix]+abs(skysubval)
  
    ;check 
    oplot, strobj.lambda[goodpix], skysubval, psym=1, color=fsc_color('blue')

    ;;write tweaked sky 
    ;savename=strmid(pixtabs,0,strpos(pixtabs,".fits"))+"_skytab.fits"
    ;mwrfits, {mastersky_wave:mastersky_wave,mastersky_flux:mastersky_flux}, savename, /crea
  
    savename=strmid(pixtabo,0,strpos(pixtabo,".fits"))+"_skycorr.fits"
    muse_writepixtab, savename, head=headobj, str=strobj
         
    ;endfor
  ;endfor
  splog, 'All done!!'

end
