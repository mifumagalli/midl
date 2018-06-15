;+
;
; Compute the slit slosses for different sources (psf, devaculeur)
;  
;
; seeing     seeing in A.S.
; slit       instrument aperture (e.g. [0.7,10]) in as
; model      structure with model parameter 
;                 1. PSF
;                 2. DeVaculeur
;
; loss       is the multiplicative factor to recover the intrinsic flux
;
;-


pro slitloss, seeing, slit, model=model, loss=loss

  
  if (model.type eq 'psf') then begin
     
     ;;prepare array
     imgps=0.01
     imgsize=1.2*max(slit)/imgps
     xarray=make_array(imgsize,/index)
     ximage=rebin(xarray,imgsize,imgsize)-floor(imgsize/2.)
     yarray=make_array(imgsize,/index)
     yimage=rebin(transpose(yarray),imgsize,imgsize)-floor(imgsize/2.)
     rimage=sqrt(ximage^2+yimage^2)


     ;;set psf model
     sigma=seeing/2.3548200/imgps 
     norm=1./(sigma^2*2*!pi)
     imgpsf=norm*exp(-(ximage^2/(2*sigma^2))-(yimage^2/(2*sigma^2)))
     
     ;;find within slit
     inside=where(ximage gt -slit[0]/2./imgps and $
                  ximage lt slit[0]/2./imgps and $
                  yimage gt -slit[1]/2./imgps and $
                  yimage lt slit[1]/2./imgps)

     ;;find slitloss
     loss=1./total(imgpsf[inside])
     
  endif
  


  if (model.type eq 'dev') then begin

     ;;set devaucouleur model at given redshift
     imgps=0.1
     
     impact_calculator, model.redshift, 1., kpc2as, /wmap7, $
                        /sil, /prop, /ang
     reffpix=model.reff*kpc2as/imgps
     
     loss=model.redshift-model.redshift
     
     ;;loop over redshift
     for zz=0, n_elements(model.redshift)-1 do begin
        

        ;;make image
        simulate_galaxy, devac, 15., reffpix[zz], 1., 0., 4., 0., 0., 0., expti=1.
        pixsize=size(devac)
        assize=pixsize[1]*imgps

        ;;trim image 
        imgsize=1.2*max(slit)/imgps
        if(imgsize lt pixsize[1]) then begin
           diff=floor((pixsize[1]-imgsize)*0.5)
           devac=devac[diff:pixsize[1]-diff,diff:pixsize[1]-diff]
        endif

        ;;convolve seeing
        devac=filter_image(devac,FWHM_GAUSSIAN=seeing/imgps)
        pixsize=size(devac)
        assize=pixsize[1]*imgps
        
        splog, "--------"
        splog, "Redshit ", model.redshift[zz]
        splog, "FOV size (as) ", assize
        splog, "Eff radius (as) ", reffpix[zz]*imgps
        ;splog, "Flux in FOV ",  -2.5*alog10(total(devac))
        splog, "--------"
    

        ;;find within slit
        xarray=make_array(pixsize[1],/index)
        ximage=rebin(xarray,pixsize[1],pixsize[1])-floor(pixsize[1]/2.)
        yarray=make_array(pixsize[1],/index)
        yimage=rebin(transpose(yarray),pixsize[1],pixsize[1])-floor(pixsize[1]/2.)
 
        inside=where(ximage gt -slit[0]/2./imgps and $
                     ximage lt slit[0]/2./imgps and $
                     yimage gt -slit[1]/2./imgps and $
                     yimage lt slit[1]/2./imgps)
        

        totflux=10^(-0.4*15.)
        insflux=total(devac[inside])
        loss[zz]=totflux/insflux
        undefine, devac
        
     endfor

  endif
  
  
end
