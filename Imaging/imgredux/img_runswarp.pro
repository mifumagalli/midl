;+
;
;  runs swarp. 
;
;-



pro img_runswarp, swarpath=swarpath, scampath=scampath, filename=filename, instr=instr,$
                  currobj=currobj, currfilt=currfilt, median=median, outname=outname, noclean=noclean, $
                  nocomb=nocomb, noresemp=noresemp


 ;;;;;;;;;;;;;;;;;;;;;;;;;;;
 ;;Run swarp to combine   ;;
 ;;;;;;;;;;;;;;;;;;;;;;;;;;;

  ;;set output subfolder
   outdir=swarpath
   ;;set location of swarp files and general stuff
   swarpconfig=getenv("MIDL")+"/Imaging/imgredux/utility/img_default.swarp"
   

   ;;create file names 
   if ~keyword_set(noresemp) then begin
      sw_img=scampath+"sci_"+filename
      sw_wgt=scampath+"wgt_"+filename
      sw_head=strarr(n_elements(filename))
 
      for ff=0, n_elements(filename)-1 do begin
         old_head=scampath+"cat_"+strmid(filename[ff],0,strlen(filename[ff])-5)+".head"
         sw_head=scampath+"sci_"+strmid(filename[ff],0,strlen(filename[ff])-5)+".head"
         ;;need to copy header
         spawn, 'cp '+old_head+" "+sw_head
      endfor
      
   endif else begin
      
      spawn, "ls "+swarpath+"/*.resamp.fits", sw_img 
      spawn, "ls "+swarpath+"/*.resamp.weight.fits", sw_wgt
      
   endelse 
   
    
   ;;get pipeline keyword
   man_general_key=['SC_OBJ','SC_GAIN','SC_RDN','SC_EXP','SC_AIR','SC_FIL','SC_SKY']
   
   ;;instrument specific
   if(instr eq 'LRIS' or instr eq 'LRISr') then begin
      man_instru_key=['DATE_BEG','DATE_END','UTC-END','TTIME','ELAPTIME','OBSERVER','TELESCOP','UTC','MJD-OBS','ST',$
                      'DATE-OBS','AIRMASS','AZ','DEC','EL','EQUINOX','HA','PONAME','RA','ROTMODE','ROTPOSN',$
                      'ROTPPOSN','REDFILT','INSTRUME','BLUFILT']
      pixel_scale = '0.135'  ;;set output pixscale - even for old ccd reproject to pixel size for blue
      size_back=128
   endif
   if(instr eq 'LBC') then begin
      man_instru_key=['EXPTIME','OBSRA','OBSDEC','OBSEPOCH','MJD_OBS','UTC_OBS','AIRMASS','INSTRUME','FILTER']
      pixel_scale = '0.224'  ;;set output pixscale
      size_back=128
   endif
   if(instr eq 'ESI') then begin
      man_instru_key=['DATE-OBS','UT','AIRMASS','RA','DEC','RAOFF','DECOFF','EQUINOX',$
                      'ROTPOSN','MJD','OBJECT','EXPOSURE','CCDGAIN','BINNING','DWFILNAM']
      pixel_scale = '0.156'  ;;set output pixscale
      size_back=128
   endif
   if(instr eq 'IMACS') then begin
      man_instru_key=['DATE-OBS','UT-TIME','AIRMASS','RA','DEC','RA-D','DEC-D','EQUINOX',$
                      'SPEED','OBJECT','EXPTIME','CCDGAIN','BINNING','EXPTYPE','FILTER','DEWAR']
      pixel_scale = '0.111'  ;;set output pixscale
      size_back=128
   endif


   
   ;;merge keyw
   ;;NOTE: swarp sees only the .head file, so all the keys have to
   ;;be set by hand.
   man_keywords=[man_general_key,man_instru_key]
   
   splog, 'Combine ', currobj, ' with filter ', currfilt
   
   ;;adjust scales from absolute to relative. Fixed to 1 now 
   newflx=fltarr(n_elements(filename))+1.

   ;;copy header info
   cphead=headfits(sw_img[0])
   man_key_val=strarr(n_elements(man_keywords))
   for kk=0, n_elements(man_keywords)-1 do $
           man_key_val[kk]=fxpar(cphead,man_keywords[kk])
   
   
   ;;create names for coadded images 
   if keyword_set(median) then begin
      outname=outdir+"sci_"+rstring(currobj)+"_"+rstring(currfilt)+"_med.fits"
      outwgt=outdir+"wgt_"+rstring(currobj)+"_"+rstring(currfilt)+"_med.fits"
   endif else begin
      outname=outdir+"sci_"+rstring(currobj)+"_"+rstring(currfilt)+".fits"
      outwgt=outdir+"wgt_"+rstring(currobj)+"_"+rstring(currfilt)+".fits"
   endelse
   
   
   ;;set swarp parameters
   if keyword_set(median) then combinetype = 'MEDIAN'  else combinetype = 'WEIGHTED'  ;;CR are with 0 weight
   if keyword_set(noclean) then cleanp='N' else cleanp='Y'
   if keyword_set(nocomb) then combp='N' else combp='Y'
   if keyword_set(noresemp) then resep='N' else resep='Y'


   splog, 'Combine using '+combinetype
   
   
   spawn, 'swarp '+strjoin(sw_img,',')+' -c '+swarpconfig+' -HEADER_SUFFIX .head -WEIGHT_IMAGE '+strjoin(sw_wgt,',')+$
          ' -IMAGEOUT_NAME '+outname+' -WEIGHTOUT_NAME '+outwgt+' -RESAMPLE_DIR '+outdir+$
          ' -GAIN_DEFAULT 1.0'+' -WRITE_FILEINFO Y'+$
          ' -WEIGHT_TYPE MAP_WEIGHT -FSCALASTRO_TYPE FIXED -FSCALE_DEFAULT '+strjoin(rstring(newflx),',')+$
          ' -SUBTRACT_BACK Y -BACK_SIZE '+rstring(size_back)+' -BACK_FILTERSIZE 3'+$
          ' -COMBINE_TYPE '+combinetype+' -CELESTIAL_TYPE NATIVE'+$
          ' -PROJECTION_TYPE TAN -PROJECTION_ERR 0.001 -CENTER_TYPE ALL -IMAGE_SIZE 0 '+$
          ' -PIXELSCALE_TYPE MANUAL -PIXEL_SCALE '+pixel_scale+' -RESAMPLING_TYPE LANCZOS3 '+$
          ' -VERBOSE_TYPE NORMAL -NTHREADS 0 -DELETE_TMPFILES '+cleanp+' -COMBINE '+combp+' -RESAMPLE '+resep, /sh
   
   ;;add the weight as extension of the fits image
   finsci=mrdfits(outname,0,h1)
   finwgt=mrdfits(outwgt,0,h2)
   
   ;;add in extensions
   mkhdr, hdr, finsci           ;generate a basic FITS header
   
   ;;add copy of keyword
   for kk=0, n_elements(man_keywords)-1 do $
           fxaddpar, hdr, man_keywords[kk], man_key_val[kk] 
   
   mwrfits, UNDEF, outname, hdr, /create
   mwrfits, finsci, outname, h1 ;;put wcs header here
   mwrfits, finwgt, outname
   
   spawn, 'rm -f '+outwgt
   
   
end
