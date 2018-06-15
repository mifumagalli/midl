;+
;
;  Given a pair of reduced pixel tables (object and sky exposures), produces
;  a SKYCORR scaled sky model for each IFU.  
;
;  pixtabo   -> name of the object pixel table (with illum corr)
;  pixtabs   -> name of the sky pixel table (with illum corr)
;  skylim   -> textfile made of 4 columns with RA and DEC limits for the sky regions
;  psky     -> percentiles of pixels used to extract the sky spectrum in the 
;              object frame  [min/max]. In the sky frame we always use [45-55].
;              if skylim is set we use [10-90] of the sky regions.
;  scalecont -> Scale the sky exposure continuum to the object exposure continuum. Useful
;               for SV data. should not be used if the sky is from a small telscope offset.
;  quick    -> select a narrow wave range to perform parameter study 
;
;-

FUNCTION SetUnion, a, b
IF a[0] LT 0 THEN RETURN, b    ;A union NULL = a
IF b[0] LT 0 THEN RETURN, a    ;B union NULL = b
RETURN, Where(Histogram([a,b], OMin = omin)) + omin ; Return combined set
END



function make_skyspec, data, lambda, nstep, npix, lambin, skysampl, wminmax, psky
    
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

       ;;loop over blocks 
       for bb=0, nbloc-2 do begin

    	  b_data=data[bb*indxbloc+offindx*steps:(bb+1)*indxbloc-1+offindx*steps]
    	  b_lambda=lambda[bb*indxbloc+offindx*steps:(bb+1)*indxbloc-1+offindx*steps]
    	  
    	  ;;find percentiles
    	  p=percentiles(b_data,value=psky)
    	  inside=where(b_data ge p[0] and b_data lt p[1],nskyp)
    	  skypixflux=[skypixflux,b_data[inside]]
    	  skypixwave=[skypixwave,b_lambda[inside]]
    			
       endfor

       splog, 'Construct spectrum for this step...'

       ;;construct a regularized version of the sky
       modelsky_flux_tmp=(modelsky_flux[*,0])*0.
       modelsky_num_tmp=(modelsky_flux[*,0])*0.
       
       skypixindex=(skypixwave-wminmax[0]+1)/skysampl
       populate_image, modelsky_flux_tmp, skypixindex, weights=skypixflux
       populate_image, modelsky_num_tmp, skypixindex
       modelsky_flux[*,steps]=modelsky_flux_tmp/modelsky_num_tmp

       ;;check 
       plot, skypixwave, skypixflux, psym=1, xrange=[6790,6890]
       oplot, modelsky_wave, modelsky_flux[*,steps], line=0, color=fsc_color('red')

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

pro muse_skytweak_v3, pixtabo, pixtabs, psky=psky, skylim=skylim, scalecont=scalecont, quick=quick, nomaster=nomaster
  
  resolve_all, /continue_on_error, /quiet
  
  ;;default 
  if ~keyword_set(psky) then psky=[0.20,0.30]
  
  ;;load pixel table 
  muse_loadpixtab, pixtabo, outhead=headobj, outstr=strobj
  muse_loadpixtab, pixtabs, outhead=headsky, outstr=strsky
  
  ;Convert to WCS positions
  aRA = fxpar(headobj.headmain, 'RA')
  aDEC = fxpar(headobj.headmain, 'DEC')
  muse_applywcs, strobj.xpix, strobj.ypix, aRA, aDEC, RApix, DECpix
  
  ;;loop over steps to fill a window
  lambin=0.4D ;;do chunks of 0.4A
  nstep=4     ;;cover the full range in 4 steps of 0.1A
  skysampl=0.4  ;;0.4 is optimal sampling in pixel table - careful with changes 
  
  if keyword_set(quick) then begin

     splog, 'Trim data...'
     
     selo=where(strobj.lambda gt 5000 and strobj.lambda lt 7000)
     sels=where(strsky.lambda gt 5000 and strsky.lambda lt 7000)
     splog, 'Selected '+strtrim(n_elements(selo),2)+' values in the pixel table.'
     
     RApix = RApix[selo]
     DECpix = DECpix[selo]
     objlambda=strobj.lambda[selo]
     objdata=strobj.data[selo]
     objdq=strobj.dq[selo]
     
     skylambda=strsky.lambda[sels]
     skydata=strsky.data[sels]
     skydq=strsky.dq[sels]
     
     undefine, selo
     undefine, sels	 
	  
  endif else begin
  
     objlambda=reform(strobj.lambda)
     objdata=reform(strobj.data)
     objdq=reform(strobj.dq)
     
     skylambda=reform(strsky.lambda)
     skydata=reform(strsky.data)
     skydq=reform(strsky.dq)
     
  endelse
  
  undefine, strsky   
  
  ;;find min max in object and sky (check if the instrument mode is consistent)
  wminmaxo=minmax(objlambda)
  wminmaxs=minmax(skylambda)
  
  if (wminmaxo[0] ne wminmaxs[0] or wminmaxo[1] ne wminmaxs[1]) then begin
   splog, 'The lambda range of the object and sky cubes does not match. Is the instrument mode the same?'
   stop
  endif else splog, 'Lambda limits: '+strtrim(wminmaxo[0],2)+'  '+strtrim(wminmaxo[1],2)
  
    ;;sort the good data to speed up loop by index 
    splog, 'Creating master object and sky spectra'
    
    if keyword_set(skylim) then begin											    
       splog, 'Reading sky limits file...'										    
       psky=[0.10,0.90] 												    
       
       readcol, skylim, RAmin, RAmax, DECmin, DECmax, comment='#', format='D,D,D,D', count=Nrect			    
    															    
       goodpixo = [-1]													    
       for rect=0, Nrect-1 do begin											    
    	 sel = where(RApix gt RAmin[rect] and RApix lt RAmax[rect] and DECpix gt DECmin[rect] and DECpix lt DECmax[rect] and objdq eq 0)  
    	 goodpixo = setunion(goodpixo, sel)											    
       endfor
       npixo = n_elements(goodpixo)													    
    endif else 	goodpixo= where(objdq eq 0 ,npixo)													    
    
    splog, 'Selected '+strtrim(npixo,2)+' values in the pixel table.'
    if npixo lt 10000 then begin
      splog, 'ERROR Too few values selected in the pixel table. Check sky masks.'
      stop
    endif					    
  
    data=temporary(objdata[goodpixo])
    lambda=temporary(objlambda[goodpixo])
  
    lsort=sort(lambda)
    data=data[lsort]
    lambda=lambda[lsort]
    
    ;;clean a bit 
    undefine, goodpixo
    undefine, lsort
    
    goodpixs=where(skydq eq 0 ,npixs)
    datas=temporary(skydata[goodpixs])
    lambdas=temporary(skylambda[goodpixs])
  
    lsort=sort(lambdas)
    datas=datas[lsort]
    lambdas=lambdas[lsort]
    
    ;;clean a bit 
    undefine, goodpixs
    undefine, lsort

    ;;make space for the sky models
    cd, current=cwd
    ifudir = cwd ;+'/IFU'+strtrim(string(ifu,format='(I02)'),2)+'_slice'+strtrim(string(slice, format='(I02)'),2)
    
    if not keyword_set(nomaster) then begin
    
       ;;run the function - needed to avoid 'edges' effect do to skylines 
       masterobj_flux =  make_skyspec( data,  lambda,  nstep, npixo, lambin, skysampl, wminmaxo, psky)
       mastersky_flux =  make_skyspec( datas, lambdas, nstep, npixs, lambin, skysampl, wminmaxo, [0.45,0.55])
       mastersky_wave =  mkarr(wminmaxo[0]-1,wminmaxo[1]+1,skysampl)
 
       ;write object and sky spectra
       savenameo=ifudir+'/mastersky_obj.txt'
       savenames=ifudir+'/mastersky_sky.txt'
       forprint, mastersky_wave, masterobj_flux, textout=savenameo, comment='#LAMBDA FLUX'
       forprint, mastersky_wave, mastersky_flux, textout=savenames, comment='#LAMBDA FLUX'
    
    endif else begin
       readcol, ifudir+'/mastersky_obj.txt', mastersky_wave, masterobj_flux
       readcol, ifudir+'/mastersky_sky.txt', mastersky_wave, mastersky_flux
    endelse
    
    ;ifuslice = strtrim(string(ifu,format='(I02)'),2)+'_'+strtrim(string(slice, format='(I02)'),2)
    ;file_mkdir, ifudir
    
    splog, 'Running skycorr...'
    
    ;prepare the SKYCORR par file and run it.
    make_skycorrpfile, ifudir, headsky.headmain
    spawn, 'skycorr skycorr.par'
    
    ;readcol, ifudir+'/mastersky_obj_SC.txt', lamcor, fluxcor, fluxSCcor, comment='#'
    ;corrsky_flux = fluxcor-fluxSCcor
    
    ;Read the results of the fit
    table = mrdfits(ifudir+'/skycorr_fit.fits',1)
    
    if keyword_set(scalecont) then begin
       splog, 'Rescale continuum spectrum...'
       wave = table.lambda
       ratiocont = table.cflux/table.mcflux

       ;clip and smooth
       okclip = where(table.class eq 0 and ratiocont gt 0.85 and ratiocont lt 1.15)

       wave = wave[okclip]
       ratiocont= median(ratiocont[okclip],10)

       ;now kill edges
       wave = wave[20:(size(wave))[1]-100]
       ratiocont = ratiocont[20:(size(wave))[1]-100]

       plot, wave, ratiocont, psym=1, yrange=[0.8,1.2]

       ;Split in two parts
       pivot=0.62
       oklow = where(wave lt pivot)
       okhigh = where(wave ge pivot)

       reslow  = poly_fit(wave[oklow], ratiocont[oklow], 3, yfit=yfit1)
       reshigh = poly_fit(wave[okhigh], ratiocont[okhigh], 3, yfit=yfit2)

       oplot, wave[oklow], yfit1, col=fsc_color('red'), thick=2
       oplot, wave[okhigh], yfit2, col=fsc_color('blue'), thick=2

       wave = table.lambda
       fitratio = fltarr(N_elements(wave))
       for i=0,n_elements(wave)-1 do begin
         if wave[i] lt pivot then fitratio[i] = reslow[0] + reslow[1] * wave[i] + reslow[2] * wave[i]^2 + reslow[3] * wave[i]^3 else $
        			  fitratio[i] = reshigh[0] + reshigh[1] * wave[i] + reshigh[2] * wave[i]^2 + reshigh[3] * wave[i]^3 
       endfor
       corrsky_flux = table.mcflux * fitratio + table.mlflux
    endif else corrsky_flux = table.mcflux  + table.mlflux
   
    splog, 'Interpolate sky over pixel table...'
    goodpix = where(strobj.dq eq 0, npixo)
    skysubval=interpol(corrsky_flux,mastersky_wave,strobj.lambda[goodpix],/spline)
  
    ;subtract
    strobj.data[goodpix]=strobj.data[goodpix]-skysubval
    strobj.stat[goodpix]=strobj.stat[goodpix]+abs(skysubval)
  
    ;check local
    plot, mastersky_wave, mastersky_flux, xrange=[6790,6890]
    oplot,mastersky_wave, masterobj_flux, color=fsc_color('red')
    oplot,mastersky_wave, corrsky_flux, color=fsc_color('blue')
    
    ;check wide
    plot, mastersky_wave, masterobj_flux-corrsky_flux, xrange=[5000,9000], yrange=[-10,10], psym=1
      
    ;Save pixel table    
    savename=strmid(pixtabo,0,strpos(pixtabo,".fits"))+"_skycorr.fits"
    muse_writepixtab, savename, head=headobj, str=strobj
         
    ;endfor
  ;endfor
  splog, 'All done!!'

end
