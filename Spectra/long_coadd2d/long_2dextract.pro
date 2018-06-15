;+
;PURPOSE
;	to cut the heart out of long_localskysub for a portable 
;	extraction program
;SYNTAX
;	extract=long_2dextract( sciimg, sciivar, waveimg, objstruct)
;
;INPUTS
;      sciimg: 	the stacked science image from long_coadd2d
;
;     sciivar: 	inverse variance
;
;     waveimg:	wavelength 2d image (with correct flexure and 
;		heliocentric correction)... output from
;		long_coadd2d
;   objstruct:	output structure from low-redux
;
;OUTPUTS
;Extract---	output structure from extraction
;  WAVE_OPT        DOUBLE    Array[NPIX]  
;	wavelength axis from optimum subtraction	
;
;  FLUX_OPT        FLOAT     Array[NPIX]
;	flux from optimum subtraction;  
;
;  SIVAR_OPT       FLOAT     Array[NPIX]
;					--
;  IVAR_OPT        FLOAT     Array[NPIX]
;					--
;  SKY_OPT         FLOAT     Array[NPIX]
;	ignore this sky tag
;					--
;  RN_OPT          FLOAT     Array[NPIX]
;  NIVAR_OPT       FLOAT     Array[NPIX]
;  MASK_OPT        BYTE      Array[NPIX]
;  FRAC_USE        FLOAT     Array[NPIX]
;  CHI2            FLOAT     Array[NPIX]
;					--
;  WAVE_BOX        DOUBLE    Array[NPIX]
;	wavelength axis from boxcar extraction
;					--
;  FLUX_BOX        FLOAT     Array[NPIX]
;	flux from boxcar extraction
;					--	
;  SIVAR_BOX       FLOAT     Array[NPIX]
;  IVAR_BOX        FLOAT     Array[NPIX]
;  NIVAR_BOX       FLOAT     Array[NPIX]
;					--
;  SKY_BOX         FLOAT     Array[NPIX]
;	ignore this sky tag
;					--
;  RN_BOX          FLOAT     Array[NPIX]
;  MASK_BOX        BYTE      Array[NPIX]
;  MINCOL          LONG          1      
;  MAXCOL          LONG          1                         
;					--
;  BOX_RAD         INT           1 
;	radius for boxcar extractions	        
;
;
;Written by R. da Silva, 8-1-09, UCSC
;-

FUNCTION long_2dextract, sciimg, sciivar, waveimg, str, prof_nsigma=prof_nsigma, box_rad=box_rad, niter=niter



if not keyword_set(niter) then niter=4L
if not keyword_set(box_rad) then box_rad=7

objstruct=str[0]
skyimage=sciimg*0


s=size(skyimage)
mask=fltarr(s[1], s[2])
xcen_med=median(objstruct.xpos)
rr=box_rad*10;  250 ;change this...
if (xcen_med - rr) LT 0 OR (xcen_med+rr) GT s[1]-1 then rr=rr/2 
mask[xcen_med-rr:xcen_med+rr, *]=1

thismask=mask
outmask=thismask


thisind = WHERE(thismask)
outmask[thisind] = ((thismask EQ 1) AND (sciivar GT 0))[thisind]


nx = (size(sciimg))[1]
ny = (size(sciimg))[2]
xarr = findgen(nx) # replicate(1.0, ny)
yarr = findgen(ny)## replicate(1.0, nx)


mincol=objstruct.mincol
maxcol=objstruct.maxcol

objwork=1
indx=[0]
modelivar=sciivar
nc = maxcol - mincol + 1L
ipix = lindgen(nc) # replicate(1, ny) + $
        replicate(1, nc) # lindgen(ny)*nx + mincol
if not keyword_set(prof_nsigma) then prof_nsigma=3.
obj_profiles=fltarr(nx*ny, objwork)

profile_struct=objstruct

 for iiter=1L, niter do begin
        splog, 'Iteration #', iiter, ' of ', niter
        img_minsky = sciimg - skyimage
        FOR ii = 0L, objwork-1L DO BEGIN
 	    iobj=0
            IF iiter EQ 1L THEN BEGIN
                splog, '-------------------REDUCING--------------'
                splog, 'Finding profile for obj #' $
                       , objstruct[indx[iobj]].OBJID, ' of ', nobj
                splog, 'At  x=', djs_median(objstruct[indx[iobj]].xpos[0]) $
                       , ' on slit #', objstruct[indx[iobj]].SLITID
                splog, 'This is Object#', indx[iobj] + 1L, ' of total ', ntot
                splog, '-----------------------------------------'
                flux = extract_boxcar(img_minsky*(mask EQ 1) $
                                      , objstruct[indx[iobj]].xpos $
                                      , objstruct[indx[iobj]].ypos $
                                      , radius = box_rad)
                mvarimg = 1.0/(modelivar + (modelivar EQ 0))
                mvar_box = extract_boxcar(mvarimg*(mask EQ 1) $
                                          , objstruct[indx[iobj]].xpos $
                                          , objstruct[indx[iobj]].ypos $
                                          , radius = box_rad)
                pixtot = extract_boxcar(float(modelivar*0 + 1) $
                                        , objstruct[indx[iobj]].xpos $
                                        , objstruct[indx[iobj]].ypos $
                                        , radius = box_rad)
                mask_box = extract_boxcar(float(mask EQ 0) $
                                          , objstruct[indx[iobj]].xpos $
                                          , objstruct[indx[iobj]].ypos $
                                          , radius = box_rad) NE pixtot
                box_denom = extract_boxcar((waveimg GT 0.0) $
                                           , objstruct[indx[iobj]].xpos $
                                           , objstruct[indx[iobj]].ypos $
                                           , radius = box_rad)
                wave = extract_boxcar(waveimg, objstruct[indx[iobj]].xpos $
                                      , objstruct[indx[iobj]].ypos $
                                     , radius = box_rad)/ $
                  (box_denom + (box_denom EQ 0))
                fluxivar = mask_box/(mvar_box + (mvar_box EQ 0))
            ENDIF ELSE BEGIN 
                last_profile = reform(obj_profiles[*, ii], nx, ny) 
                trace = objstruct[indx[iobj]].xpos ## replicate(1, nx)
                objmask =  ((xarr GE (trace - 2.0*box_rad)) AND $
                            (xarr LE (trace + 2.0*box_rad)))
                temp_str = long_extract_optimal(thismask*waveimg $
                                                , img_minsky $
                                                , sciivar*thismask $
                                                , last_profile $
                                                , outmask*objmask  $
                                                , skyimage $
                                                , objstruct[indx[iobj]].xpos $
                                                , modelivar =  $
                                                modelivar*thismask $
                                                , box_rad = box_rad $
                                                , rn_img = rn_img)
                ;; Bad extraction?
                if temp_str.wave_opt[0] GE 0. then begin
                    ;; Good extraction
                    flux = temp_str.FLUX_OPT
                    fluxivar = temp_str.IVAR_OPT
                    wave = temp_str.WAVE_OPT
                endif
             ENDELSE

            obj_profiles[ipix, ii] = long_gprofile(img_minsky[ipix] $
                                                   , (modelivar*mask)[ipix] $
                                                   , waveimg[ipix] $         
                                                   , objstruct[indx[iobj]].xpos $
                                                   - mincol $
                                                   , wave, flux, fluxivar $
                                                   , objstruct[indx[iobj]] $
                                                   , hwidth = $
                                                   objstruct[indx[iobj]].maskwidth $
                                                   , fwhmfit = fwhmfit $
                                                   , thisfwhm = $
                                                   objstruct[indx[iobj]].fwhm $
                                                   , xnew = xnew, nccd = nccd $
                                                   , SN_GAUSS = SN_GAUSS $
                                                   , PROF_NSIGMA= $
                                                   prof_nsigma[indx[iobj]] $
                                                   , MED_SN2 = MED_SN2, /silent)
            ;; Hand apertures are currently not included in the final
            ;; object model image
;            ignoreobj[ii] = (objstruct[indx[iobj]].HAND_AP EQ 1)
;objstruct[indx[iobj]].xpos    = xnew + mincol
            objstruct[indx[iobj]].fwhmfit = fwhmfit
            objstruct[indx[iobj]].fwhm    = median(fwhmfit)
            IF iiter EQ niter THEN BEGIN
               ;; We are not writing these out to disk at present
	        
                profile_struct={OBJID:objstruct[indx[iobj]].objid,$	
                SLITID:objstruct[indx[iobj]].slitid, $
                MED_SN2:MED_SN2, $ 
		fwhm:median(fwhmfit), $
		fwhmfit:fwhmfit,$
                PROFILE:reform(obj_profiles[*, ii], nx, ny)}
	        temp_str=struct_addtags(profile_struct, temp_str)		
;                IF KEYWORD_SET(profile_filename) THEN $                
;                  mwrfits, profile_struct, profile_filename $
;                                  , CREATE = (indx[iobj] EQ 0), /silent
            ENDIF
        ENDFOR
ENDFOR

return, temp_str

end
