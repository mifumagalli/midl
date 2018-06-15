;#############################################################################
;
; MF- > the main driver for ppxf: adapted from example of Cap. 
;
; Usage example for the procedure PPXF, which
; implements the Penalized Pixel-Fitting (pPXF) method by
; Cappellari M., & Emsellem E., 2004, PASP, 116, 138.
; The example also shows how to include a library of templates
; and how to mask gas emission lines if present.
; 
;
;  gal_lin, noi_lin, wave  input spectrum in linear units. Wave has to
;                          be rest frame and equally spaced 
; 
;  FWHM_gal                spectral FWHM in rest frame 
;  fit*                    output paramters erros in the fit 
;  bias                    passed on to fit  
;

pro ppxf_kinematics_drive, gal_lin, noi_lin, wave, FWHM_gal=FWHM_gal, fitsolution=fitsolution,$
                           fiterror=fiterror, bias=bias 
  

  ;;define spectrum edges (MF)
  lamRange1 = minmax(wave)
  
  ;;inputs are assumed to be in rest frame already 
  ;;
  ;;z = 0.01544                   ; Initial estimate of the galaxy redshift
  ;;lamRange1 = lamRange1/(1+z)   ; Compute approximate restframe wavelength range
  ;;FWHM_gal = FWHM_gal/(1+z)     ; Adjust resolution in Angstrom


  ;;normalize the spectrum in linear space 
  noi_lin = noi_lin/median(gal_lin)
  gal_lin = gal_lin/median(gal_lin) 
  
  ;;rebin in log 
  ppxf_log_rebin, lamRange1, gal_lin, galaxy, logLam1, VELSCALE=velScale
  ppxf_log_rebin, lamRange1, noi_lin^2, noise, logLam1err, VELSCALE=velScale
  noise=sqrt(noise)


  ;; Read the list of filenames from the Single Stellar Population library
  ;; by CENARRO 
  cenarro = file_search(getenv('MIDL')+'/IFU/voronppxf/cenarro_models/scan*.fits',COUNT=nfiles)
  FWHM_tem = 1.5                ; Vazdekis spectra have a resolution FWHM of 1.8A.

  ;; Extract the wavelength range and logarithmically rebin one spectrum
  ;; to the same velocity scale of the input galaxy spectrum, to determine
  ;; the size needed for the array which will contain the template spectra.
  
  fits_read, cenarro[0], ssp, h2
  lamRange2 = sxpar(h2,'CRVAL1') + [0d,sxpar(h2,'CDELT1')*(sxpar(h2,'NAXIS1')-1d)]
  ppxf_log_rebin, lamRange2, ssp, sspNew, logLam2, VELSCALE=velScale
  templates = dblarr(n_elements(sspNew),nfiles)
  
  ;; Quadratic sigma difference in pixels Cenarro --> instrument 
  ;; The formula below is rigorously valid if the shapes of the 
  ;; instrumental spectral profiles are well approximated by Gaussians. 
  FWHM_dif = SQRT(FWHM_gal^2 - FWHM_tem^2)
  sigma = FWHM_dif/2.355/sxpar(h2,'CDELT1') ; Sigma difference in pixels 

  ;; IMPORTANT: To avoid spurious velocity offsets of the templates, the
  ;; NPIXEL keyword in PSF_GAUSSIAN must be an odd integer as done below
  
  lsf = psf_Gaussian(NPIXEL=2*ceil(4*sigma)+1, ST_DEV=sigma, /NORM, NDIM=1)
  for j=0,nfiles-1 do begin
     fits_read, cenarro[j], ssp
     ssp = convol(ssp,lsf)      ; Degrade template to instrument resolution
     ppxf_log_rebin, lamRange2, ssp, sspNew, VELSCALE=velScale
     templates[*,j] = sspNew/median(sspNew) ; Normalizes templates 
  endfor
  
  
; The galaxy and the template spectra do not have the same starting wavelength.
; For this reason an extra velocity shift DV has to be applied to the template
; to fit the galaxy spectrum. We remove this artificial shift by using the
; keyword VSYST in the call to PPXF below, so that all velocities are
; measured with respect to DV. This assume the redshift is negligible.
; In the case of a high-redshift galaxy one should de-redshift its 
; wavelength to the rest frame before using the line below (see above).
;
  c = 299792.458d
  dv = (logLam2[0]-logLam1[0])*c ; km/s
  
  vel = 1d                      ; Initial estimate of the galaxy velocity in km/s
  goodPixels = ppxf_determine_goodPixels(logLam1,lamRange2,vel)
  
  ;; Here the actual fit starts. The best fit is plotted on the screen.
  ;;Gas emission lines are excluded from the pPXF fit using the GOODPIXELS keyword.
  
  start = [vel, 180d]           ; (km/s), starting guess for [V,sigma]

  ppxf, templates, galaxy, noise, velScale, start, sol, GOODPIXELS=goodPixels, $
        /PLOT, MOMENTS=4, DEGREE=4, VSYST=dv, ERROR=error, bias=bias 

  splog, 'Formal errors:    dV    dsigma       dh3       dh4'
  splog, error[0:3]*sqrt(sol[6]), FORMAT='(10x,2f10.1,2f10.3)'

  ;;pass back results [Vel,Sigma,h3,h4,h5,h6,Chi^2/DOF]
  fitsolution=sol[0:3]
  fiterror=error[0:3]*sqrt(sol[6])


end
;------------------------------------------------------------------------------
