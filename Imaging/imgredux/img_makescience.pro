;+
;procedure to make science frame (one each time)
;bias,zero,dark,flat,gain,exptime
;
;name    --> single image to process
;filter  --> filter of the image
;STATUS  --> the status structre from ccdproc 
;gzip    --> compress final science frame
;path    --> data path
;object  --> object name to update header
;instr   --> the instrument used
;time    --> the exposure time of the image
;
;
;-

pro img_makescience, name, filter, time, status=status,$
                     skyfit=skyfit, gzip=gzip, root=root, $
                     side=side, path=path, object=object, instr=instr, notweakgain=notweakgain, nocr=nocr
  
  common ccdinfo, xfull,yfull,nchip,namp,biaslev,satur, goodata
  
  splog, "Processing sci frame ", name
  head=headfits(path+strtrim(name,2),exten=0)
 
  
  ;;open file (instrument dependent) 
  if(instr eq 'LRIS') then begin
     date=fxpar(head,"DATE")
     if(date gt '2009-04-01') then chips=lris_readfits(path+strtrim(name,2),header=header,/notrim,layout=layout)
     if(date lt '2009-04-01') then chips=lris_readold(path+strtrim(name,2),header=header,layout=layout)
  endif
  if(instr eq 'LRISr') then chips=lris_readold(path+strtrim(name,2),header=header,layout=layout)
  if(instr eq 'LBC') then lbtc_readfits,path+strtrim(name,2),mosaic=chips,/notrim,layout=layout,/nodisp,header=header
  if(instr eq 'ESI') then begin
     chips=esi_readfits(path+strtrim(name,2),header=header,layout=layout)
     date=fxpar(head,"DATE")
  endif
  if(instr eq 'FEINC') then chips=feinc_readfits(path+strtrim(name,2),header=header,layout=layout)
  if(instr eq 'IMACS') then chips=imacs_readfits(path+strtrim(name,2),header=header,layout=layout)

  
  ;;go for oscan subtraction
  img_oscan, chips, imageout, layout=layout, instr=instr
  
  ;;remove zero if set
  if keyword_set(bias) then imageout=imageout-bias
  
  ;;find corresponding flat field and flat
  flat=where(strcompress(status.filflat,/remove_all) eq strcompress(filter,/remove_all))
  nameflat=side+root+"_"+strcompress(status.filflat[flat],/remove_all)+"flat.fits"
  flatfits=mrdfits('proc/'+nameflat[0],/silent)     
  
  ;;apply flat
  divideflat, imageout, flatfits, invvar=invvar, minval=0.01, /quiet
    
  ;;apply gain (this is instrument dependet)
  if(instr eq 'LRIS') then img_lris_gain, imageout, gimage, outgain=gain_value, $
          outscale=outscale, side=side,notweakgain=notweakgain, date=date
  if(instr eq 'LRISr') then img_lrisr_gain, imageout, gimage, outgain=gain_value, $
          outscale=outscale, side=side
  if(instr eq 'LBC') then img_lbtc_gain, imageout, gimage, outgain=gain_value, $
          outscale=outscale, side=side
  if(instr eq 'ESI') then img_esi_gain, imageout, gimage, header, outgain=gain_value, $
          outscale=outscale
  if(instr eq 'FEINC') then img_feinc_gain, imageout, gimage, header, outgain=gain_value, $
          outscale=outscale
  if(instr eq 'IMACS') then img_imacs_gain, imageout, gimage, header, outgain=gain_value, $
          outscale=outscale, side=side
  
  ;;ivar image (gain has been already applied to science)
  img_makeweight, gimage, flatfits, weight, time=time, $
                  instr=instr, side=side, readnoise=readnoise, date=date, name=name
  
  ;;flag cosmic rays using lacosmic (if it works, ok, if not too bad...)
  ;;requires to write to disk 2 images
  
  ;;tweak the flags on specific instrument
  ;;no difference for now in each side        
  if(instr eq 'LRIS' and side eq 'R') then thenlac_sigma=4. else thenlac_sigma=5. 
  lac_iter=3.
 
  if(~keyword_set(nocr)) then begin
     if(instr ne "FEINC") then begin ;;skip cr for some instrument
        
        splog, 'Flag CR...'
        
        ;;make random suffix
        randid=string(randomu(seed,1)*1000.,format='(I04)')
        
        ;;write to disk
        mwrfits, gimage, 'proc/'+"tmp_crimg"+randid+".fits", header, /CREATE, /SILENT
        la_cosmic, ['proc/'+"tmp_crimg"+randid+".fits"], sigclip=lac_sigma, gain=1., niter=lac_iter, $
                   masksuffix="-mask", outsuff="-out", readn=readnoise
        
        ;;read in mask and make a good mask
        inmask=mrdfits('proc/'+"tmp_crimg"+randid+"-mask.fits",0,/sil)
        flag=where(inmask gt 0,nflg)
        spawn, "rm -f "+'proc/'+"tmp_crimg"+randid+".fits "+'proc/'+"tmp_crimg"+randid+"-mask.fits"
        ;;add CR to weight mask
        if(nflg gt 0) then weight[flag]=0.
        
     endif
  endif

  ;;per second
  gimage=gimage/time
  weight=weight*time^2
  
  ;;update header 
  sxaddpar, header, "OBJECT", object 
  sxaddpar, header, "RED_DATE", string("Redux ended ", SYSTIME())  
  sxaddpar, header, "COMMENT", "IMG reduced using IMG_CCDPROC", AFTER='RED_DATE'
  sxaddpar, header, "COMMENT", string("Flat field ",nameflat), AFTER='RED_DATE'
  sxaddpar, header, "COMMENT", string("Image in electron/second. Diveded by T=",time), AFTER='RED_DATE'
  sxaddpar, header, "COMMENT", strjoin(string("Gain applied: ",gain_value)), AFTER='RED_DATE'
  sxaddpar, header, "COMMENT", strjoin(string("Scale applied: ",outscale)), AFTER='RED_DATE'
  
  ;;for lbt update bzero
  if(instr eq 'LBC' or instr eq 'IMACS') then sxaddpar, header, "BZERO", 0.0
  if(instr eq 'IMACS') then begin
    sxdelpar, header, "BIASSEC"
    sxdelpar, header, "DATASEC"
  endif

  ;;write final image, weight and mask
  savename=strmid(name,0,strpos(name,'.fits'))+"_redux.fits"
  mwrfits, gimage, 'proc/'+savename, header, /CREATE, /SILENT
  mwrfits, weight, 'proc/'+savename
  
  if keyword_set(gzip) then begin
     ;;gzip image
     command="gzip "+'proc/'+savename
     spawn, command
  endif else splog, "Done with "+savename


end
