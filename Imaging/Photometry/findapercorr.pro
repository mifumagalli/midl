;+
;
;  This is an interactive tools that allows you to 
;  normalize the flux level for the psf model and extracts 
;  a value for aperture correction at 1.35FWHM with relative errors.
;  All the fwhm units are normalized.
;
;  model       the txt file from modelpsf.pro
;  corr        in output, the aperture correction at 2FWHM
;  err         in output, an estimate of the normalization error
;              for Delta_fwhm=1
;  fwhmnorm    in output, the fwhm used for normalisation
;  fluxnorm    in output, the flux used for normalisation
;
;
;-


pro findapercorr, model, corr=corr, err=err, fwhmnorm=fwhmnorm, $
                  fluxnorm=fluxnorm



;;load the file
readcol, model, Rad, Flux, Derivative, /sil


quit='y'
repeat begin
    ;;plot the flux
    window, 0
    plot, rad, flux, psym=1, xtitle='Radius (FWHM units)', ytitle='Relative Flux enclosed'
    
    ;;plot the derivative
    window, 1
    plot, rad, derivative, psym=1, xtitle='Radius (FWHM units)', ytitle='Error/FWHM', /ylog
    
    ;;extract the position
    cursor, fwhmnorm, y
    splog, 'Normalize at ', fwhmnorm
    wdelete, 0, 1
    

    
    ;;noramalize and find aperture correction
    fluxnorm=interpol(flux,rad,fwhmnorm)
    array_flux_norm=flux/fluxnorm
    corr=1./interpol(array_flux_norm,rad,1.35)
    err=interpol(derivative,rad,fwhmnorm)

    splog, 'Correction ', corr 
    splog, 'Error ', err

    ;;prompt whether it's good 
    window, 0
    plot, rad, array_flux_norm, psym=1, xtitle='Radius (FWHM units)', ytitle='Relative Flux enclosed', $
      yrange=[0,1.1]
    oplot, [1.35], [1./corr], psym=4
    oplot, [-100000,100000], [1.,1.],  line=1
    oplot, [1.35,1.35],[-100000,100000], line=2

    
    read, 'Is this ok? [y/n]', quit
    
endrep until quit eq 'y'











end
