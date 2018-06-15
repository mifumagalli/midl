;+
;
;
; For an image and a source structure (even with multi filters)
; find good stars and create a psf model 
;
; fits      --> the data array
; objstr    --> structure of the sources    
; seg       --> array of the segmentation map
; filt      --> for multifilter stucture, the index of the filter 
; stack     --> set to a name where the fits file of the stacked psf
;               is saved plus a gaussian 2D model
; model     --> if set, a bspline elliptical model and residual to fits file. 
; plot      --> set to a filename to save the ps of the PSF details
; nstar     --> nstar to use (too few will produce a non
;               representative median, too many will bring in galaxies)
; mask      --> any mask /wght image that can be used to mask defects (0=reject)
; mean      --> if set, do a mean stack. Defualt is median
; txt       --> if set to name the flux enclosed is saved in a file
; manstar   --> if set, one can inspect and skip stars:
;               y: include in stack 
;               n: do not include in stack
;               s: skip and generate stack
; sdss      --> if there are sdss info in the structure, use sdss stars
;-



pro modelpsf, fits, objstr, seg, filt=filt, model=model, stack=stack, plot=plot, $
              nstar=nstar, mask=mask, mean=mean, txt=txt, sdss=sdss

;;default
if ~keyword_set(filt) then filt=0
if ~keyword_set(stack) then stack='psf_stack.fits'
if ~keyword_set(plot) then plot='psf_model.ps'
if ~keyword_set(txt) then txt='psf_model.txt'

if ~keyword_set(nstar) then nstar=25

;;find good stars getting rid of defects
if keyword_set(sdss) then $
  nodef=where(objstr.defect[filt] lt 1 and objstr.flags[filt] lt 4 and $
              objstr.fwhm_world[filt]*3600 lt 4*objstr.seeing[filt] and $
              objstr.sdss_star[filt] gt 0) else $
  nodef=where(objstr.defect[filt] lt 1 and objstr.flags[filt] lt 4 and $
              objstr.fwhm_world[filt]*3600 lt 4*objstr.seeing[filt]) 


star=reverse(sort(objstr[nodef].class_star[filt]))
select=star[0:nstar-1]

;;select box radius
maxfwhm=max(objstr[nodef[select]].fwhm_image[filt])
radius=8.2*maxfwhm


for i=0, nstar-1 do begin 
;;select center
xpos=objstr[nodef[select[i]]].xwin_image[filt]
ypos=objstr[nodef[select[i]]].ywin_image[filt]

;;extract box 
fits_box=extractbox(fits,xpos,ypos,radius)
seg_box=extractbox(seg,xpos,ypos,radius)
;xatv, fits_box, /bl

;;take care of the mask
if keyword_set(mask) then begin
    mas_box=extractbox(mask,xpos,ypos,radius)
    bad=where(mas_box eq 0,nbad) 
    if(nbad gt 0) then fits_box[bad]=0
endif

;;additional sky sub
noobj=where(seg_box lt 1)
mmm, fits_box[noobj], skyval, skysigma
fits_box=temporary(fits_box-skyval)

;;flag other objects (but keep sky)
objid=objstr[nodef[select[i]]].number[0] ;;obj id stored only in 1st position
mas=where((seg_box gt 0 and seg_box lt objid) or (seg_box gt objid),nm)
if(nm gt 0)then fits_box[mas]=0.0
;xatv, seg_box, /bl

;;stack normalized
fits_box=temporary(fits_box/total(fits_box))
;xatv, fits_box, /bl


if keyword_set(mean) then begin 
    ;;mean
    if(i eq 0.) then stackf=fits_box else stackf=stackf+fits_box
endif else begin
    ;;median
    if(i eq 0.) then begin
        ;;set stack
        st_size=size(fits_box)
        stackf=fltarr(nstar,st_size[1],st_size[2])
    endif 
    stackf[i,0:st_size[1]-1,0:st_size[2]-1]=fits_box
endelse

endfor

if ~keyword_set(mean) then stackf=djs_median(stackf,1)


;;derive some statistics used for photometry
;;with gaussian fit
gmodel=gauss2dfit(stackf,gfit,/tilt) 

splog, ' Residual sky ', gfit[0]
splog, ' Normalization ', gfit[1]
splog, ' sigma_x ', gfit[2]
splog, ' sigma_y ', gfit[3]
splog, ' PA CCW from X (deg) ', gfit[6]*180./!PI

;;find approximation for 2D FWHM and 1D FWHM
fwhm2d=2.*sqrt(2.*alog(2.))*sqrt(gfit[2]^2+gfit[3]^2) 
fwhm1d=2.*sqrt(2.*alog(2.))*0.5*(gfit[2]+gfit[3]) 

splog, ' 2D FWHM (approx) ', fwhm2d 
splog, ' 1D FWHM (approx) ', fwhm1d


;;save fits file with gaussian model
mwrfits, stackf, stack, /create
mwrfits, gmodel, stack
size=size(stackf)

;;find light enclosed in fraction of 1d_fwhm
;;The normalisation is very sensitive to the chosen sky value
core=mkarr(0.1,2.,0.2)*fwhm1d
wing=mkarr(3.,8.,0.7)*fwhm1d

aper, stackf, gfit[4], gfit[5],  flux_core, errap, sky, skyerr, 1., core, [1,2], $
  [1,1], /EXACT, /FLUX, /SILENT, SETSKYVAL=0.
aper, stackf, gfit[4], gfit[5],  flux_wing, errap, sky, skyerr, 1., wing, [1,2], $
  [1,1], /EXACT, /FLUX, /SILENT, SETSKYVAL=0.


;;set the normalisation arbitrary. Later on one can change it 
;;inspecting the flux derivative 
flux=[flux_core,flux_wing]
rad_nfw=[core,wing]/fwhm1d
derivative=(flux-shift(flux,1))/(rad_nfw-shift(rad_nfw,1))
derivative[0]=(flux[1]-flux[0])/(rad_nfw[1]-rad_nfw[0])

;;normalize
flux_n=flux/flux[n_elements(flux)-1]

    ;;plot and print
m_psopen, plot, /land
plot, rad_nfw, flux_n, psym=1, title=' 1D FWHM (pix) '+string(fwhm1d,format='(F5.2)'),$
  xtitle='Radius (FWHM units)', ytitle='Relative Flux enclosed'
plot, rad_nfw, derivative, psym=1, title=' 1D FWHM '+string(fwhm1d,format='(F5.2)'),$
  xtitle='Radius (FWHM units)', ytitle='Error/FWHM', /ylog
m_psclose

forprint, rad_nfw, flux, derivative, textout=txt, comment='Rad (FWHM),Flux(absolute),Derivative'

;;make a bspline model
if keyword_set(model) then begin
;;;This needs more work!!
;     size=size(stackf)
;     ;;position image    
;     xx = reverse(findgen(size[1])) # replicate(1.,size[1]) - float(0.5*size[1])
;     yy = replicate(1.,size[2]) # findgen(size[2]) - float(0.5*size[2])
;     delta_xx= xx - gfit[4]
;     delta_yy= yy - gfit[5]
;     r=sqrt(delta_xx^2+delta_yy^2)
;     theta=atan(delta_yy,delta_xx)
;     ;;fake invvar
;     invvar=0.*stackf-stackf+0.01
;     ;;multipole parameters
;     ntheta=[0,-1,1,-2,2,-4,4]      ; monop, dipol, etc... for fit with psf convolution
;     rbkpt=findgen(10)              ; breakpoints every pix within 10 FWHM
;     rset=bspline_radial(r, theta, stackf, invvar=invvar, $
;                         rbkpt=rbkpt, ntheta=ntheta, yfit=bmodel,upper=-100, lower=100) 
;     ;bmod=bspline_ellip(par, x=x, y=y, data=data, invvar=invvar, $
;     ;                   psf=psf, deviates=deviates, ntheta=ntheta, rbkpt=rbkpt, sset=sset)
endif


end
