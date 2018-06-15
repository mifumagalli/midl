;+
;
;  Given a fully processed cube with no sky subtraction, construct a
;  local sky  
;
;  pixtab   -> name of the pixel table for skysub (with illum corr)
;  cubename -> if set to name of the cube, apply sky and produce
;              images for test of sky level  
;  psky     -> percentiles of pixels used to model the sky  [min/max]
;  quick    -> select a narrow wave range to perform parameter study 
;
;-

pro muse_skymodel, pixtab, psky=psky, cubename=cubename, quick=quick

  ;;default 
  if ~keyword_set(psky) then psky=[0.20,0.30]
  
  ;;load pixel table 
  muse_loadpixtab, pixtab, outhead=head, outstr=str
  
  ;;loop over steps to fill a window
  lambin=0.4D ;;do chunks of 0.4A
  nstep=4     ;;cover the full range in 4 steps of 0.1A


  if keyword_set(quick) then begin

     splog, 'Trim data...'
     sel=where(str.lambda gt 6200 and str.lambda lt 6900)
     str.lambda=str.lambda[sel]
     str.data=str.data[sel]
     str.dq=str.dq[sel]
     str.stat=str.stat[sel]
     
     str={lambda:str.lambda[sel],data:str.data[sel],dq:str.dq[sel],stat:str.stat[sel],$
          xpix:str.xpix[sel],ypix:str.ypix[sel],origin:str.origin[sel]}
     
  endif
  
  ;;find min max
  wminmax=minmax(str.lambda)
  
  ;;construct wave array of sky model
  skysampl=0.4  ;;0.4 is optimal sampling in pixel table - careful with changes 
  modelsky_wave=mkarr(wminmax[0]-1,wminmax[1]+1,skysampl)
  
  ;;sort the good data to speed up loop by index 
  splog, 'Sort good data... '
  goodpix=where(str.dq eq 0,npix)
  data=str.data[goodpix]
  var=str.stat[goodpix]
  lambda=str.lambda[goodpix]
  
  lsort=sort(lambda)
  data=data[lsort]
  var=var[lsort]
  lambda=lambda[lsort]

  ;;clean a bit 
  undefine, lsort
  
  ;;find blocks of wave - do it by index which is way faster
  nbloc=floor((wminmax[1]-wminmax[0])/lambin)
  indxbloc=floor(npix/nbloc)
  offindx=1.*indxbloc/nstep

  ;;make space for the sky models
  modelsky_flux=fltarr(n_elements(modelsky_wave),nstep)
  
  ;;run the loop - needed to avoid 'edges' effect do to skylines 
  for steps=0, nstep-1 do begin
     
     splog, 'Iteration ', steps+1, ' out of ',  nstep
     
     ;;create storage for pixels
     skypixflux=0.
     skypixwave=0.
     skypixvar=0.

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
        skypixvar=[skypixvar,b_var[inside]]
                      
     endfor

     splog, 'Construct sky model for this step...'

     ;;construct a regularized version of the sky
     modelsky_flux_tmp=modelsky_wave*0.
     modelsky_num_tmp=modelsky_wave*0.
     
     skypixindex=(skypixwave-wminmax[0]+1)/skysampl
     populate_image, modelsky_flux_tmp, skypixindex, weights=skypixflux
     populate_image, modelsky_num_tmp, skypixindex
     modelsky_flux[*,steps]=modelsky_flux_tmp/modelsky_num_tmp

     ;;check 
     plot, skypixwave, skypixflux, psym=1, xrange=[6790,6890]
     oplot, modelsky_wave, modelsky_flux[*,steps], line=1, color=fsc_color('red')
     
     ;;do some checks on cube 
     if keyword_set(cubename) then begin
        
        splog, 'Run checks on cube'

        ;;load cube -- use obs wave at this stage
        muse_loadcube, cubename, wave=cwave, flux=cflux, var=cvar, obswave=cobswave

        ;select wave 
        sz=size(cflux)
        cwsel=where(cobswave ge wminmax[0] and cobswave le wminmax[1],ncwsel)

        ;;eval sky 
        cskyval=interpol(modelsky_flux[*,steps],modelsky_wave,cobswave[cwsel],/spline)
        
        ;;apply sky to cube in wave range
        for cww=0, ncwsel-1 do begin
           cflux[*,*,cwsel[cww]]=cflux[*,*,cwsel[cww]]-cskyval[cww]
        endfor

        ;;collapse cube over wrange
        white=median(cflux[*,*,cwsel],dime=3)
        ;xatv, white, /bl
        
        ;;unwrap cube 
        ifu_unwrap, cobswave[cwsel],(mk_finite(cflux[*,*,cwsel])), twodspc
        ;xatv, twodspc, /bl

        stop

     endif
     

  endfor
  
  splog, 'Construct master sky model...'
  
  mastersky_flux=median(modelsky_flux,dim=2)
  mastersky_wave=modelsky_wave
 
  ;;check 
  plot, mastersky_wave, mastersky_flux, xrange=[6790,6890]
  for steps=0, nstep-1 do oplot, modelsky_wave, modelsky_flux[*,steps], line=1, color=fsc_color('red')
  
  splog, 'Interpolate sky over pixel table...'
  skysubval=interpol(mastersky_flux,mastersky_wave,str.lambda[goodpix],/spline)
  
  ;;subtract
  ;str.data[goodpix]=str.data[goodpix]-skysubval
  ;str.stat[goodpix]=str.stat[goodpix]+abs(skysubval)
  
  ;;check 
  ;;oplot, str.lambda[goodpix], skysubval, psym=1, color=fsc_color('blue')
  
  splog, 'Write output...'

  ;;write sky 
  savename=strmid(pixtab,0,strpos(pixtab,".fits"))+"_skytab.fits"
  mwrfits, {mastersky_wave:mastersky_wave,mastersky_flux:mastersky_flux}, savename, /crea
   
  ;;write pixel table 
  ;savename=strmid(pixtab,0,strpos(pixtab,".fits"))+"_sky.fits"
  ;muse_writepixtab, savename, head=head, str=str
  
  splog, 'All done!!'

end
