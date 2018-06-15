;+
;
;   Run a combine only version of swarp.
;
;
;   suffix -> name use to create file list 
;   headref  -> file name with path used to copy some basic head keyowrd
;
;-


pro img_mancomb, suffix, instr=instr, headref=headref
  
  
  if ~keyword_set(instr) then instr="LRIS"
  
  spawn, "ls *"+suffix, list
  
  
  ;;set location of swarp files and general stuff
  swarpconfig=getenv("MIDL")+"/Imaging/imgredux/utility/img_default.swarp"
  
  ;;set swarp parameters
  weightopts = ' -WEIGHT_TYPE MAP_WEIGHT -WEIGHT_SUFFIX .weight.fits'
  fluxscaleopts = ' -FSCALE_KEYWORD SW_SCALE -FSCALASTRO_TYPE FIXED' ;;some error at the edge of field
  combinetype = 'WEIGHTED' ;;CR are with 0 weight
  coordtype = 'EQUATORIAL'
  projection_err = '0.001'
  centeropts = ' -CENTER_TYPE ALL -IMAGE_SIZE 0'
  
  ;;get pipeline keyword
  man_general_key=['SC_OBJ','SC_GAIN','SC_RDN','SC_EXP','SC_AIR','SC_FIL','SC_SKY']
  
  ;;instrument specific
  if(instr eq 'LRIS' or instr eq 'LRISr') then begin
      man_instru_key=['DATE_BEG','DATE_END','UTC-END','TTIME','ELAPTIME','OBSERVER','TELESCOP','UTC','MJD-OBS','ST',$
                      'DATE-OBS','AIRMASS','AZ','DEC','EL','EQUINOX','HA','PONAME','RA','ROTMODE','ROTPOSN',$
                      'ROTPPOSN','REDFILT','INSTRUME','BLUFILT']
  endif
  if(instr eq 'LBC') then begin
      man_instru_key=['EXPTIME','OBSRA','OBSDEC','OBSEPOCH','MJD_OBS','UTC_OBS','AIRMASS','INSTRUME','FILTER']
  endif
  if(instr eq 'ESI') then begin
      man_instru_key=['DATE-OBS','UT','AIRMASS','RA','DEC','RAOFF','DECOFF','EQUINOX',$
                      'ROTPOSN','MJD','OBJECT','EXPOSURE','CCDGAIN','BINNING','DWFILNAM']
  endif
  
  ;;merge keyw
  man_keywords=[man_general_key,man_instru_key]
  
  ;;NOTE: swarp sees only the .head file, so all the keys have to
  ;;be set by hand.
  
  outname="combined.fits"
  outwgt="combined.wgt.fits"
  
  ;;copy header info
  cphead=headfits(headref)
  man_key_val=strarr(n_elements(man_keywords))
  for kk=0, n_elements(man_keywords)-1 do $
          man_key_val[kk]=fxpar(cphead,man_keywords[kk])
  

  spawn, 'swarp '+strjoin(list,',')+' -c '+swarpconfig+$
          ' -IMAGEOUT_NAME '+outname+' -WEIGHTOUT_NAME '+outwgt+$
          weightopts+' -HEADER_SUFFIX .head'+fluxscaleopts+$
          ' -COMBINE_TYPE '+combinetype+' -CELESTIAL_TYPE '+coordtype+$
          ' -PROJECTION_TYPE TAN -PROJECTION_ERR '+projection_err+centeropts+$
          ' -GAIN_KEYWORD SC_GAIN -GAIN_DEFAULT 1.0'+$
          ' -COMBINE_BUFSIZE 1024 -WRITE_FILEINFO Y'+$
          ' -VERBOSE_TYPE NORMAL -NTHREADS 4 -DELETE_TMPFILES Y -COMBINE Y -RESAMPLE N', /sh
  
  
  ;;add the weight as extension of the fits image
  finsci=mrdfits(outname,0,h1)
  finwgt=mrdfits(outwgt,0,h2)
  
  ;;add in extensions
  mkhdr, hdr, finsci            ;generate a basic FITS header
 

  ;;add copy of keyword
  for kk=0, n_elements(man_keywords)-1 do $
          fxaddpar, hdr, man_keywords[kk], man_key_val[kk] 
  
  mwrfits, UNDEF, outname, hdr, /create
  mwrfits, finsci, outname, h1 ;;put wcs header here
  mwrfits, finwgt, outname
  
  spawn, 'rm -f '+outwgt
  spawn, 'rm -f *.resamp.*'

  
end
