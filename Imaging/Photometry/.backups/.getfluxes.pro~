;+
;
; This procedure take a science image, a segmenation map  
; from sexctractor, an object structure, a noise model,  
; and compute aperture photometry in kron radii. Colors are  
; computed using the segmentation map. For smaller object, 
; circular apertures and aperture correction is used.
;
; scimap    the fits file of the science assumed in e-/s
; segmap    the segmenation map
; objstr    object structure with sext parameters
; weight    extract using a weight/mask map to mask out bad pixels
; noise     a noise model [sigma1,alpha,beta]. For uncorrelated noise
;           use alpha=1,beta=1. sigma=sigma_1*alpha*n_pix^beta 
; texp      the effective exposure time (used for error estimate)
; zp        photometric zero point
; fwhm_mod  fwhm in pixel from the model. Photometry for objects which are smaller
;           than 2*FWHM is computed in circular aperture with 1.35fwhm
;           aperture in diameter (optimal extraction).
; apercorr  if provided, an aperture correction is applied.        
; errapcor  if provided, this 
; filt      or multifilter structure, this is the index of the filter 
;
; fluxout     output   
; fluxerrout  output
; localsky    output
; final_id    output
;
;
;+


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;Subroutine that does the full photometry
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

pro gf_dophot, index, filt, objstr, fit_box, fwhm_mod, $
               apercorr, errapcor, noise, texp,  kron=kron

;;set aperture parameters
if keyword_set(kron) then begin
    pos_ang=objstr[index].theta_image[filt]+90
    a_ell=objstr[index].kron_radius[filt]*objstr[index].a_image[filt]
    b_ell=objstr[index].kron_radius[filt]*objstr[index].b_image[filt]
endif else begin
    pos_ang=0.
    a_ell=1.35*fwhm_mod
    b_ell=1.35*fwhm_mod
endelse

;;make image
bxsiz=size(fit_box)
bxpos=0.5*floor(size(fit_box))
excen=bxpos[1]
eycen=bxpos[2]


;;;;;;;;;The box is centered on the window position. Avoid recentering
;;cntrd, fit_box, bxpos[1], bxpos[2], excen, eycen,  fwhm_mod, /KEEPCENTER, /SILENT
;;;;fix if fails
;;excen=(mk_finite(excen))[0] 
;;eycen=(mk_finite(eycen))[0] 
;;if(excen lt 0) then  excen=bxpos[1]
;;if(eycen lt 0) then  eycen=bxpos[2]


dist_ellipse, ellimage, [bxsiz[1],bxsiz[2]], excen, eycen, a_ell/b_ell, pos_ang
fluxpix=where(ellimage le a_ell, nfpx)
if(nfpx gt 0) then fluxtot=total(fit_box[fluxpix]) else fluxtot=0.
        
;;compute noise with positive flux    
posflux=fluxtot > 0.
err_fluxtot=sqrt((noise[0]*noise[1]*nfpx^noise[2])^2+posflux/texp)

;;store
objstr[index].tot_flux[filt]=fluxtot
objstr[index].errtot_flux[filt]=err_fluxtot

if keyword_set(kron) then begin 
;;do not correct
    objstr[index].tot_flux_cor[filt]=fluxtot
    objstr[index].errtot_flux_cor[filt]=err_fluxtot
    ;;set SN
    objstr[index].tot_sn[filt]=posflux/err_fluxtot

endif else begin
    ;;correct
    objstr[index].tot_flux_cor[filt]=apercorr*fluxtot
    objstr[index].errtot_flux_cor[filt]=sqrt((apercorr*err_fluxtot)^2+errapcor^2)
    ;;set SN
    objstr[index].tot_sn[filt]=apercorr*posflux/sqrt((apercorr*err_fluxtot)^2+errapcor^2)
    
endelse


;print, posflux, objstr[index].flux_auto[filt]
;xatv,  fit_box,  /blo
;xatv,  ellimage,  /blo
  
end


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;Subroutine that extract the box and 
;does the mask and sky sub
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
pro gf_extract, index, filt, objstr, fwhm_mod, fit_box, seg_box,$
                noise, scimap, segmap, weight, idthis, skyval, color=color, $
                contin=contin

;;find object location 
;;For color, you should use (x_ and y_ coordinates) to preserve
;;the aligment. For photometry, use x_win and y_win
if keyword_set(color) then begin
    xpos=objstr[index].x_image[filt]
    ypos=objstr[index].y_image[filt]
endif else begin
    xpos=objstr[index].xwin_image[filt]
    ypos=objstr[index].ywin_image[filt]
endelse

;;cut a box aroudn the object for sky evaluation and analysis
;;Take 2 times the equivalent radius of isophotal area 
;;or 10 times the fwhm
boxrad=2*sqrt(objstr[index].isoarea_image[filt]) > 10.*fwhm_mod

;;deal with much smaller images to speed this up
fit_box=extractbox(scimap,xpos,ypos,boxrad)
seg_box=extractbox(segmap,xpos,ypos,boxrad)
wgh_box=extractbox(weight,xpos,ypos,boxrad)

;;turn the wgt into mask
mas=where(wgh_box gt 0,nmk)
mas_box=wgh_box
if(nmk gt 0) then mas_box[mas]=1.

;;mask other objects but keep sky
idthis=objstr[index].number[0]
mas=where((seg_box gt 0 and seg_box lt idthis) or (seg_box gt idthis),nm)
if(nm gt 0) then mas_box[mas]=0.

;;if the object is too close to the edge, there may be no good
;;pixels at all. In this case, stop working on this...
;;At least 20 pixel are needed for sky

if(total(mas_box) gt 25) then begin
    
    contin=1
    
    ;;evaluate local sky excluding defects and objects
    ;;If there are too few pixels, mmm gives an error and 
    ;;the code crashes. You should force mmm by adding 
    ;;/continue and a return after the error message
    skyind=where(mas_box gt 0 and seg_box lt 1,nsk)
    if(nsk gt 25) then mmm, fit_box[skyind], skyval, skysigma else begin
        skyval=0.
        skysigma=0.
        ;;if no good pixels flag as non-go object
        if(nsk le 0) then contin=0
    endelse
    
    objstr[index].sky_local[filt]=skyval
    fit_box=temporary(fit_box-skyval)
    
    ;;apply mask
    fit_box=fit_box*mas_box
    
    ;;if noise is not set, find noise information from sky
    if(noise[0] lt 0) then noise[0]=skysigma
    
endif else begin
    contin=0
    skyval=0
endelse

;free
undefine, wgh_box, mas_box

;;xatv,  fit_box,  /blo

end


pro getfluxes,  scimap, segmap, objstr, noise=noise, texp=texp, weight=weight, $
                zp=zp, fluxout=final_flux, fluxerrout=final_fluxerr, $
                localsky=localsky, final_id=final_id, fwhm_mod=fwhm_mod, apercorr=apercorr,$
                errapcor=errapcor, filt=filt

   
;;ii=reverse(sort(objstr.flux_iso[filt]))
;;objstr=objstr[ii]


;;set some default
if ~keyword_set(noise) then noise=[-1.,1.,1.] ;uncorrelated noise
if ~keyword_set(texp) then texp=1.
if ~keyword_set(zp) then zp=0.
if ~keyword_set(filt) then filt=0.
if ~keyword_set(fwhm_mod) then fwhm_mod=5.
if ~keyword_set(apercorr) then apercorr=1.
if ~keyword_set(errapcor) then errapcor=0.

;;get image info
size=size(scimap)

;;make equal good pixel image if not provided
if ~keyword_set(weight) then weight=scimap-scimap+1.

;;get the number of object which are not marked as defects by 
;;sextractor and make space for output
indobj=where(objstr.flags[filt] lt 4,nobj)
defect=where(objstr.flags[filt] ge 4,ndef)

final_mag=fltarr(nobj)
final_magerr=fltarr(nobj)
localsky=fltarr(nobj)
final_id=fltarr(nobj)

;;get kron parameter
minrad=min(objstr[where(objstr.kron_radius[filt] gt 0)].kron_radius[filt])

splog, 'Compute mag for ', nobj, ' objects'
splog, 'Start at ', systime()

;;loop over all the objects
for oj=0L, nobj-1 do begin

    
    ;;cut a first box for color determination
    gf_extract, indobj[oj], filt, objstr, fwhm_mod, fit_box, seg_box,$
      noise, scimap, segmap, weight, idthis, skyval, contin=contin, /color

    ;;store id
    final_id[oj]=idthis
 
    
    if(contin gt 0) then begin
        
        ;;color 
        objloc=where(seg_box eq idthis,npixcol)
        ;;compute color flux and error
        if(npixcol gt 0) then colorflux=total(fit_box[objloc]) else colorflux=0.
        ;;if flux is negative use only back error, which is the dominant error
        posflux=colorflux > 0.
        colorfluxerr=sqrt((noise[0]*noise[1]*npixcol^noise[2])^2+posflux/texp)
        ;;store
        objstr[indobj[oj]].color_flux[filt]=colorflux
        objstr[indobj[oj]].errcolor_flux[filt]=colorfluxerr
        ;;set detection signal to noise
        objstr[indobj[oj]].color_sn[filt]=posflux/colorfluxerr
 
        ;;free 
        undefine, fit_box, seg_box
    endif
    
    ;;now cut a second box for total flux determination
    ;;Here gets set the sky level in the final structure
    gf_extract, indobj[oj], filt, objstr, fwhm_mod, fit_box, seg_box,$
      noise, scimap, segmap, weight, idthis, skyval, contin=contin
    
    ;;store sky
    localsky[oj]=skyval

    if(contin gt 0) then begin
    
        ;;total flux
        if(objstr[indobj[oj]].kron_radius[filt] ge minrad and $
           objstr[indobj[oj]].fwhm_image[filt] gt 1.*fwhm_mod) then begin
            
            ;;extract with kron
            gf_dophot, indobj[oj], filt, objstr, fit_box, fwhm_mod, $
              apercorr, errapcor, noise, texp, /kron
            
        endif else  gf_dophot, indobj[oj], filt, objstr, fit_box,$
          fwhm_mod, apercorr, errapcor, noise, texp ;;circular aperture
        
        ;;free
        undefine, fit_box, seg_box
    endif 

    
;;;;;;;;;;;;;;;;OLD
;set ellipse parameters 
;     xc=objstr[oj].x_image
;     yc=objstr[oj].y_image
;     kron aperture   NOT SURE IF THIS IS THE CORRECT KRON APERTURE!!!!!
;     pos_ang=objstr[oj].theta_image+90
;     a_ell=objstr[oj].kron_radius*objstr[oj].a_image
;     b_ell=objstr[oj].kron_radius*objstr[oj].b_image
;    
;     limit to min radius (take into account kron_radius=0 for
;     flux<0)
;     if(a_ell lt minrad*objstr[oj].a_image or b_ell lt minrad*objstr[oj].b_image) then begin
;         a_ell=minrad*objstr[oj].a_image
;         b_ell=minrad*objstr[oj].b_image
;     endif
;    
;     dist_ellipse, ellimage, [xsize,ysize], xc, yc, a_ell/b_ell, pos_ang
;;;;;;;;;;;;;;;;;;;;
     
endfor

splog, 'End at ', systime()


end
