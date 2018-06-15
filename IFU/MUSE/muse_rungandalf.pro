;
; Run Gandalf without dereddening option to subtract the stellar 
; continuum spectrum from each spaxel
; 
; cubefile = name of the input datacube
; image = image of the FoV to be used to define where to perform the subtraction
; limits = an array like [x1,x2,y1,y2,fl1,fl2] where fl is the flux threshold in the image
;          for X and Y both the extremes are included
; outfile = the name of the output file
; silent = limit the verbosity to the progress messages. 
;          As a default all the gandalf output is printed
; ebv_gal = the value of E(B-V) for the Galaxy
;

;-----------------------------------------------------------
pro my_statusline, step_count, abs_step, t_start

step_count += 1

if step_count mod abs_step eq 0 then begin
   curr_percent = 5*step_count/abs_step
   t_here = systime(/seconds)
   if curr_percent lt 100 then begin
    eta = ((t_here-t_start)/(curr_percent)*(100-curr_percent))
    case 1 of
      eta le 60 : etastr=strtrim(fix(eta),2)+'s'
      eta gt 60 and eta le 3600 : etastr=strtrim(fix(eta/60),2)+':'+strtrim(string(fix(eta mod 60),format='(I02)'),2)+' mm:ss'
      eta gt 3600 : etastr=strtrim(fix(eta/3600),2)+':'+strtrim(string(fix((eta mod 3600)/60),format='(I02)'),2)+':'+$
                    strtrim(string(fix((eta mod 3600) mod 60),format='(I02)'),2)+' hh:mm:ss'  
    endcase
    print, '[PROGRESS] '+strtrim(curr_percent,2)+'% processed - ETA: '+etastr 
   endif else begin
    tot = (t_here-t_start)
    case 1 of
      tot le 60 : totstr=strtrim(fix(tot),2)+'s'
      tot gt 60 and tot le 3600 : totstr=strtrim(fix(tot/60),2)+':'+strtrim(string(fix(tot mod 60),format='(I02)'),2)+' mm:ss'
      tot gt 3600 : totstr=strtrim(fix(tot/3600),2)+':'+strtrim(string(fix((tot mod 3600)/60),format='(I02)'),2)+':'+$
                    strtrim(string(fix((tot mod 3600) mod 60),format='(I02)'),2)+' hh:mm:ss'  
    endcase
    print, '[PROGRESS] '+strtrim(curr_percent,2)+'% processed - Total time: '+totstr 
   endelse 
endif

end

;----------------------------------------------------------------------
function mask_emission_lines,npix,Vsys,emission_setup,velscale,l0_gal,lstep_gal,$
                             sigma=sigma,l_rf_range=l_rf_range,log10=log10, silent=silent

; Return a list of goodpixels to fit that excludes regions potentially
; affected by gas emission and by sky lines. Unless the log10 keyword
; is specified, wavelength values are assumed to be ln-rebinned, and
; are defined by the l0_gal, lstep_gal, npix parameters. The position of
; gas and sky emission lines is set by the input emission_setup
; structure and the width of the mask by the sigma parameter. If a
; sigma value is not passed than the width of each line is taken from
; the emission_setup structure.
;
; The rest-frame fitting wavelength range can be manually restricted
; using the l_rf_range keyword to pass min and max observed
; wavelength. Typically used to exclude regions at either side of
; spectra.


; speed of light
c = 299792.458d
; define good pixels array
goodpixels = range(0,npix-1) 
; if set, exclude regions at either ends of the spectra using the keyword l_rf_range
if keyword_set(l_rf_range) then begin
    pix0     = ceil((alog(l_rf_range[0])-l0_gal)/lstep_gal+Vsys/velscale)
    pix1     = ceil((alog(l_rf_range[1])-l0_gal)/lstep_gal+Vsys/velscale)
    if keyword_set(log10) then begin
        pix0     = ceil((alog10(l_rf_range[0])-l0_gal)/lstep_gal+Vsys/velscale)
        pix1     = ceil((alog10(l_rf_range[1])-l0_gal)/lstep_gal+Vsys/velscale)
    endif
    goodpixels = range(max([pix0,0]),min([pix1,npix-1])) ; NEW - V1.3 
endif

tmppixels  = goodpixels
; looping over the listed emission-lines and mask those tagged with an
; 'm' for mask. Mask sky lines at rest-frame wavelength
for i = 0,n_elements(emission_setup.i)-1 do begin
    if (emission_setup.action[i] eq 'm') then begin
        if NOT KEYWORD_SET(silent) then print,'--> masking ' + emission_setup.name[i]
        if (emission_setup.name[i] ne 'sky') then $
          meml_cpix = ceil((alog(emission_setup.lambda[i])-l0_gal)/lstep_gal+Vsys/velscale)
        if (emission_setup.name[i] ne 'sky' and keyword_set(log10)) then $
          meml_cpix = ceil((alog10(emission_setup.lambda[i])-l0_gal)/lstep_gal+Vsys/velscale)
        ; sky lines are at rest-frame
        if (emission_setup.name[i] eq 'sky') then $
          meml_cpix = ceil((alog(emission_setup.lambda[i])-l0_gal)/lstep_gal) 
        if (emission_setup.name[i] eq 'sky' and keyword_set(log10)) then $
          meml_cpix = ceil((alog10(emission_setup.lambda[i])-l0_gal)/lstep_gal) 
        ; set the width of the mask in pixels using either
        ; 3 times the sigma of each line in the emission-line setup 
        ; or the provided sigma value 
        if keyword_set(sigma) then msigma = 3*sigma/velscale
        if not keyword_set(sigma) then msigma = 3*emission_setup.s[i]/velscale
        meml_bpix = meml_cpix - msigma
        meml_rpix = meml_cpix + msigma
        w = where(goodpixels ge meml_bpix and goodpixels le meml_rpix) 
        if (w[0] ne -1) then begin
            tmppixels[w] = -1 
        endif else begin 
            if NOT KEYWORD_SET(silent) then print,'this line is outside your wavelength range. We shall ignore it' 
            emission_setup.action[i] = 'i'                                         ; NEW - V1.3 
        endelse
    endif
endfor
w = where(tmppixels ne -1)
goodpixels = goodpixels[w]

return,goodpixels
end

function remouve_detected_emission,galaxy,bestfit,emission_templates,sol_gas_A,AoN_thresholds,$
                                   AoN=AoN,goodpixel=goodpixel
; Given the galaxy spectrum, the best fit, the emission-line
; amplitudes, a vector with the A/N threshold for each line, and the
; array containing the spectra of each best-fitting emission-line
; templates, this function simply compute the residual-noise lvl,
; compares it the amplitude of each lines, and remouves from the
; galaxy spectrum only the best-matching emission-line templates for
; which the correspoding A/N exceed the input threshold.  This is a
; necessary step prior to measurements of the strength of the stellar
; absorption line features
;
; A list of goodpixels may be optionally input, for instance if any
; pixel was excluded by sigma-clipping during the continuum and
; emission-line fitting, or by excluding pixels on either ends of the
; spectra
;
; Also optionally outputs the computed A/N ratios.

; Get the Residual Noise lvl.
resid = galaxy - bestfit
if keyword_set(goodpixel) then resid = resid[goodpixels]
resid_noise = robust_sigma(resid, /ZERO)
; A/N of each line
AoN = sol_gas_A/resid_noise
; Create neat spectrum, that is, a spectrum where only detected
; emission has been remouved
neat_galaxy = galaxy
for i = 0,n_elements(sol_gas_A)-1 do begin
    if (AoN[i] ge AoN_thresholds[i]) then neat_galaxy = neat_galaxy - emission_templates[*,i]
endfor

return,neat_galaxy
end

function dust_calzetti,l0_gal,lstep_gal,npix,ebv,vstar,LOG10=log10
; This procedure uses the dust model of Calzetti et al. (2000, ApJ,
; 533, 682), and for a given E(B-V) value returns the flux attenuation
; array, which can be used to get reddened templates. Here the spectra
; are assumed to be binned on a ln-rebinned wavelentgh grid as defined
; by input l0_gal,lstep_gal,npix parameters. The input receiding
; velocity vstar, is used to derive the dust reddening in the galaxy
; rest-frame.
; 
; Can be used also to de-reddened the galaxy spectra by the Milky-Way
; dust extinction, using as E(B-V) the opposite of the Schlegel et
; al. values found in NED and vstar = 0.
;
; Initial version kindly provided by S. Kaviray, Oxford, 2006.

; reconstruct the wavelength array in Anstroms, and compute rest-frame
; values
lambda = exp(dindgen(npix)*lstep_gal + l0_gal)
if keyword_set(log10) then lambda = 10^(dindgen(npix)*lstep_gal + l0_gal)
lambda = lambda/exp(vstar/299792.458d)

; array to hold k(lambda) values
k = fltarr(n_elements(lambda))           

for i=0,n_elements(lambda)-1 do begin
     ; convert wavelength units from angstroms to micrometres
     l = lambda(i)/1e4                   
     ; assign k values
     if (l ge 0.63 and l le 2.2) then k(i) = 2.659*(-1.857+1.040/l)+4.05
     if (l lt 0.63)              then k(i) = 2.659*(-2.156+1.509/l-0.198/l^2+0.011/l^3)+4.05
     if (l gt 2.2)               then k(i) = 0.0
endfor

return,(10^(-0.4*ebv*k))
; this should be then multiplied by the spectrum flux array
end

;-----------------------------------------------------
;MAIN PROGFRAM
pro muse_gandalf, cubefile, image, limits, outfile, plot=plot, silent=silent, ebv_gal=ebv_gal

resolve_all, /quiet, /continue_on_error

;read datacube and variance cube
print,'Reading the datacube ...'

datacube   = mrdfits(cubefile, 1, datahead)
statcube   = mrdfits(cubefile, 2, errhead)

image = mrdfits(image,0,headima)


; Reading the templates
; 
; These must have a) the same log-lambda step EQUAL to that of the
; galaxy spectrum AND b) the same spectral resolution of the galaxy
; spectrum. 

print,'Reading template library ...'

vazdekis = file_search('miles_models/Mun1.30Z*.fits',COUNT=nfiles)
FWHM_tem = 2.51 ; Vazdekis+10 spectra have a resolution FWHM of 2.51A.

; Extract the wavelength range and logarithmically rebin one data and one template spectrum
lamrange    = sxpar(datahead, 'CRVAL3') + [0,sxpar(datahead, 'CD3_3')*sxpar(datahead, 'NAXIS3')]
log10_rebin, lamrange, reform(datacube[100,100,*]), galaxy, log10lam, velscale=velocity

fits_read, vazdekis[0], ssp, hdr_templ
lamRange2 = sxpar(hdr_templ,'CRVAL1') + [0d,sxpar(hdr_templ,'CDELT1')*(sxpar(hdr_templ,'NAXIS1')-1d)]
log10_rebin, lamRange2, ssp, sspNew, logLam2, VELSCALE=velocity
n_templ=nfiles
n_pix_templ = n_elements(sspNew)
templates = dblarr(n_pix_templ,n_templ)

print,'Template array created!'
FOR j=0, nfiles-1 DO BEGIN
    fits_read, vazdekis[j], ssp
    log10_rebin, lamRange2, ssp, sspNew, VELSCALE=velocity
    templates[*,j] = sspNew/median(sspNew) ; Normalizes templates 
END

;Prepare other variables
redshift = 0.01555
wavesol_lin = sxpar(datahead, 'CRVAL3') + sxpar(datahead, 'CD3_3')*findgen(sxpar(datahead, 'NAXIS3'))
wavesol_log = 10^(alog10(sxpar(datahead, 'CRVAL3'))-0.000045 + ((log10lam[1] - log10lam[0])*findgen(sxpar(datahead, 'NAXIS3'))))
okHalin = where(wavesol_lin gt 6600 and wavesol_lin lt 6900)
okHalog = where(wavesol_log gt 6600 and wavesol_log lt 6900)
okHblin = where(wavesol_lin gt 4800 and wavesol_lin lt 5100)
okHblog = where(wavesol_log gt 4800 and wavesol_log lt 5100)
plotpos = cgLayout([1,2], OXMargin=[3, 3], OYMargin=[3, 3], XGap=5, YGap=5)

c = 299792.4580d ; Speed of light in km/s
mdegree=3
t_start    = systime(/seconds)
abs_step = ((long(1+limits[1]-limits[0])*(1+limits[3]-limits[2]) * 5)/100)
step_count = 0L

print,'Starting the big loop...go and grab a coffee and some sandiwches'
for col=limits[0],limits[1] do begin
  for row=limits[2],limits[3] do begin
    
    my_statusline, step_count, abs_step, t_start
    if image[col,row] gt limits[4] and image[col,row] lt limits[5] then begin
       
       spec = reform(datacube[col,row,*])
       err  = reform(statcube[col,row,*])
       
       log10_rebin, lamrange, spec, galaxy, log10lam, velscale=velocity; ,/flux
       log10_rebin, lamrange, err, error, log10lam, velscale=velocity;   ,/flux

       ; which need to be adjusted to exclude zeros
       error = sqrt(error)
       w = where(error eq 0) & if (w[0] ne -1) then error[w] = 1e6

       ;check rebin
       ;plot, wavesol_lin[okHalin], spec[okHalin], col=fsc_color('white')
       ;oplot, wavesol_log[okHalog], specnew[okHalog], col=fsc_color('red')

       ;GANDALF_SDSS,INFILE='temp_MUSE_rebin.fits', $
       ;  EMISSION_SETUP_FILE='emission_lines_setup_with_Balmer_decrement',$
       ;  MDEGREE=3 ;,/PLOT , EBV_GAL=0.20
       

       IF KEYWORD_SET(plot) THEN window,0,xsize=500,ysize=500, retain=2

       ; Set the intrinsic velocity dispersion resolution of the data.
       int_disp = 49.5

       ; B) De-redden the spectra for the galactic extinction in the
       ; direction of the target, using Schlegel et al values from NED
       IF KEYWORD_SET(ebv_gal) THEN BEGIN
           l0_gal   = salog10(sxpar(datahead, 'CRVAL3'))-0.000045
           lstep_gal = log10lam[1] - log10lam[0]
           dereddening_attenuation = DUST_CALZETTI(l0_gal,lstep_gal,n_elements(galaxy),-ebv_gal,0.0d,/log10)
           galaxy = galaxy*dereddening_attenuation
       ENDIF

       ; E) Computing the velocity offset between the starting wavelengths
       ; of the galaxy and template spectra. You need to include a alog(10)
       ; factor since we are using log10 rebinned data instead of ln rebinned.
       if NOT KEYWORD_SET(silent) then print,'--> set template-galaxy spectra offset'
       l0_gal	=  alog10(sxpar(datahead, 'CRVAL3'))-0.000045
       l0_templ =  alog10(sxpar(hdr_templ,'CRVAL1'))
       offset = -(l0_gal-l0_templ)*c*alog(10.0d)


       ; F) Preamble to PPXF
       if NOT KEYWORD_SET(silent) then print,'--> and initial guesses for pPXF using the SDSS redshifts'
       ; this is the velocity scale spanned by each log_lambda interval
       ; You need to include a alog(10) factor as well.
       lstep_gal = log10lam[1] - log10lam[0]
       velscale  = c*lstep_gal*alog(10.0d)
       ; Initial V and sigma guesses, in km/s. For V we use the receiding
       ; velocity derived by the SDSS (in the low-Z regime !!!), plus the
       ; previously derived velocity offset. For sigma we stick at least to
       ; 150 km/s
       SDSS_z = redshift
       SDSS_s = -99
       if (SDSS_s gt velscale/2.0d) then $
         start  = [alog(SDSS_z+1)*c+offset,SDSS_s] else $
         start  = [alog(SDSS_z+1)*c+offset,150.0d]

       ; G) masking emission-lines
       ; read in emission-line setup files
       if NOT KEYWORD_SET(silent) then print,'--> reading in the emission-line setup file, and creating the corresponding structure'
       eml_file = 'emission_lines_setup_with_Balmer_decrement'
       readcol,eml_file,eml_i,eml_name,eml_lambda,eml_action,eml_kind,eml_a,eml_v,eml_s,eml_fit,$
         f='(i,a,f,a,a,f,f,f,a)',skipline=2,comment='#',/silent
       ; creating emission setup structure, to be used for masking the
       ; emission-line contaminated regions and fitting the emission lines in GANDALF
       emission_setup = create_struct('i',eml_i,'name',eml_name,'lambda',eml_lambda,'action',eml_action,$
        			      'kind',eml_kind,'a',eml_a,'v',eml_v,'s',eml_s,'fit',eml_fit)

       ; call the function mask_emission_lines, which loops over the listed
       ; emission-lines and mask those tagged with an 'm' for mask. Use
       ; emission-line masks with a sigma of 150 km/s, and use only the
       ; wavelength region corresponding to the rest-frame extent of the
       ; templates.
       l_rf_range = 10^[l0_templ,l0_templ+n_pix_templ*(log10lam[1] - log10lam[0])]		       ; NEW - V1.3
       goodpixels = mask_emission_lines(n_elements(galaxy),alog(SDSS_z+1)*c,emission_setup,velscale,$  ; NEW - V1.3
        				l0_gal,lstep_gal,sigma=150,/log10,l_rf_range=l_rf_range, /silent)        ; NEW - V1.3

       ; H) PPXF fit! Fit only V and sigma
       ; A constant additive polynomial (degree=0) is used in together with
       ; the multiplicative polynomials (always recommended).
       PPXF, templates, galaxy, error, velscale, start, sol, bestfit=bestfitppxf, $
           goodpixels=goodpixels,plot=plot, moments=2, degree=0, mdegree=mdegree, quiet=silent
       if NOT KEYWORD_SET(silent) then  print,'--> Fit to the stellar continuum masking regions potentially affected by gas emission!'
       if keyword_set(DEBUG) then pause

       ; J) Preamble to GANDALF
       ; Switch the tag 'action' from 'm' (for mask) to 'f' (for fit) for all
       ; lines we wish to fit, i.e, do not touch the ones with an
       ; 'i'. Likewise keep masking the regions affected by interestella NaD
       ; absorption and which were originally affected by sky lines.
       emission_setup_orig = emission_setup
       i_lines = where(emission_setup.action eq 'm' and emission_setup.name ne 'sky') ; NEW - V1.3
       emission_setup.action[i_lines] = 'f'

       ; Re-assign the goodpixels array, masking what is still tagged with an 'm'
       goodpixels = mask_emission_lines(n_elements(galaxy),alog(SDSS_z+1)*c,emission_setup,velscale,$
        				l0_gal,lstep_gal,sigma=150,/log10,l_rf_range=l_rf_range, /silent) ; NEW - V1.3
       if NOT KEYWORD_SET(silent) then print,'--> lift the emission-line mask...'

       ; Prepare emission_setup structure for GANDALF, which should only
       ; deal with the lines we fit
       i_f = where(emission_setup.action eq 'f') 
       dummy = emission_setup
       emission_setup = create_struct('i',dummy.i[i_f],'name',dummy.name[i_f],$
        			      'lambda',dummy.lambda[i_f],'action',dummy.action[i_f],$
        			      'kind',dummy.kind[i_f],'a',dummy.a[i_f],$
        			      'v',dummy.v[i_f],'s',dummy.s[i_f],$
        			      'fit',dummy.fit[i_f])
       ; Assign the stellar systemic velocity as initial guess for the gas
       ; kinematics, adding also the velocity offset specified in the
       ; emission-line setup. This makes it easier to add blue or red wings to
       ; the line profile.
       emission_setup.v = sol[0]-offset + emission_setup.v ; NEW - V1.3
       ; save the stellar kinematics
       sol_star = sol

       ; H) Call Gandalf, giving it only the stellar kinematics as input
       ; sol. Now include reddening 

       if NOT KEYWORD_SET(silent) then print,'--> and fit simultaneously the stellar continuum and the ionised-gas emission lines!'
       GANDALF, templates, galaxy, error, velscale, sol, emission_setup, $
         l0_gal, lstep_gal, GOODPIXELS=goodpixels, INT_DISP=int_disp, $
         BESTFIT=bestfit, EMISSION_TEMPLATES=emission_templates, WEIGHTS=weights, $
         PLOT=plot, /LOG10, L0_TEMPL=l0_templ, quiet=silent, $
         FOR_ERRORS=for_errors, ERROR=esol, mdegree=mdegree ;, REDDENING=[0.05,0.05]

       if keyword_set(DEBUG) then pause

       ; K) Make the unconvolved optimal template. This can be useful for the
       ; line-strength analysis in order to work out the necessary correction
       ; due measured indices due to kinematic broadening.
       ;nweights = weights[0:n_templ-1]/total(weights[0:n_templ-1])
       ;otemplate = dblarr(n_elements(templates[*,0]))
       ;for j=0,n_templ-1 do otemplate = otemplate + templates[*,j]*nweights[j]
       ;print,'--> Creating the unconvolved optimal template spectrum...'


       ; L) Call the routine that remouves only the detected emission,
       ; assuming a constant A/N cut of 4. The routine will also output the
       ; A/N of each line.

       ; Best fitting amplitudes
       ;i_l = where(emission_setup.kind eq 'l') 
       ;sol_gas_A = sol[dindgen(n_elements(i_l))*4+1]
       ; constant A/N=4 threshold
       ;AoN_thresholds = dblarr(n_elements(i_l)) + 4.0
       ;spec_neat = remouve_detected_emission(galaxy,bestfit,emission_templates,sol_gas_A,AoN_thresholds,AoN=sol_gas_AoN)
       ;print,'--> ... and cleaning the galaxy spectrum from any detected gas emission line'

       ; M) Add to the emission setup structure the flux, amplitude and
       ; kinematics of each line, and call the result the fit_results
       ; structure. Add to it also the A/N values, the stellar kinematics,
       ; and the normalised template weights.  This is to save not only the
       ; emission-line fitting results but also the conditions under which
       ; the fit was performed.
       ;i_l = where(emission_setup.kind eq 'l') 
       ;sol_gas_F = sol[dindgen(n_elements(i_l))*4+0]
       ;sol_gas_A = sol[dindgen(n_elements(i_l))*4+1]
       ;sol_gas_V = sol[dindgen(n_elements(i_l))*4+2]
       ;sol_gas_S = sol[dindgen(n_elements(i_l))*4+3]
       ;
       ;sol_EBmV  = sol[n_elements(i_l)*4:*]

       ;dummy = emission_setup
       ;if keyword_set(FOR_ERRORS) then begin
       ;    esol_gas_F = esol[dindgen(n_elements(i_l))*4+0]
       ;    esol_gas_A = esol[dindgen(n_elements(i_l))*4+1]
       ;    esol_gas_V = esol[dindgen(n_elements(i_l))*4+2]
       ;    esol_gas_S = esol[dindgen(n_elements(i_l))*4+3]
       ;    
       ;    esol_EBmV  = esol[n_elements(i_l)*4:*]
       ;
       ;    fit_results = create_struct('i',dummy.i[i_l],'name',dummy.name[i_l],'lambda',dummy.lambda[i_l],$
       ; 			       'action',dummy.action[i_l],'kind',dummy.kind[i_l],'a',dummy.a[i_l],$
       ; 			       'v',dummy.v[i_l],'s',dummy.s[i_l],'fit',dummy.fit[i_l],$
       ; 			       'Flux',sol_gas_F,'Ampl',sol_gas_A,'Vel',sol_gas_V,'Sigma',sol_gas_S,$
       ; 			       'eFlux',esol_gas_F,'eAmpl',esol_gas_A,'eVel',esol_gas_V,'eSigma',esol_gas_S,$
       ; 			       'AoN',sol_gas_AoN,'EBmV',sol_EBmV,'eEBmV',esol_EBmV,'Vel_stars',sol_star[0]-offset,$
       ; 			       'Sigma_stars',sol_star[1],'Norm_Weights',nweights)
       ;endif else begin
       ;    fit_results = create_struct('i',dummy.i[i_l],'name',dummy.name[i_l],'lambda',dummy.lambda[i_l],$
       ; 			       'action',dummy.action[i_l],'kind',dummy.kind[i_l],'a',dummy.a[i_l],$
       ; 			       'v',dummy.v[i_l],'s',dummy.s[i_l],'fit',dummy.fit[i_l],$
       ; 			       'Flux',sol_gas_F,'Ampl',sol_gas_A,'Vel',sol_gas_V,'Sigma',sol_gas_S,$
       ; 			       'AoN',sol_gas_AoN,'EBmV',sol_EBmV,'Vel_stars',sol_star[0]-offset,$
       ; 			       'Sigma_stars',sol_star[1],'Norm_Weights',nweights)
       ;endelse


       ;print,'--> Finally, write out the fit results!'
       ; N) Finally writing out the outputs

       ; writing fits file with the fitting results structure
       ;pars_file = (strsplit(infile,'.',/extract))[0]+'_pars.fits'
       ;mwrfits,fit_results,pars_file,/create,/silent

       ; writing multi-extension fits file with fit results
       ;fits_file = (strsplit(infile,'.',/extract))[0]+'_fits.fits'
       ;writefits,fits_file,galaxy    ,hdr_gal
       ;writefits,fits_file,bestfit   ,/append ; appending best fitting model
       
       ; Add up the best-fitting emission templates to get the emission spectrum
       if ((size(emission_templates))[0] eq 1) then emission = emission_templates	   
       if ((size(emission_templates))[0] eq 2) then emission = total(emission_templates,2)
       
       ;writefits,fits_file,emission  ,/append ; appending best-fitting emission-line spectrum
       ;writefits,fits_file,spec_neat ,/append ; appending galaxy spectrum without detected emission
       ;writefits,fits_file,bestfitppxf  ,/append ; appending best fitting model

       ; NEW - V1.3
       ;dummy = galaxy*0 & dummy[goodpixels] = 1 & goodpixels = dummy
       ;writefits,fits_file,goodpixels ,/append ; appending goodpixels array (1 good, 0, masked)
       ; writing fits file with optimal template
       ;otempl_file = (strsplit(infile,'.',/extract))[0]+'_otempl.fits'
       ;writefits,otempl_file,otemplate ,hdr_templ 
       
       best_fit = bestfit
       line_fit = emission
       cont_fit = best_fit-line_fit
       
       if NOT KEYWORD_SET(silent) then print, 'Interpolate model over the spectrum wave sol...'
       cont_fit_lin=interpol(cont_fit,wavesol_log,wavesol_lin,/spline)
      
       ;check plot
       plot, wavesol_lin[okHalin], spec[okHalin], col=fsc_color('white'), position=plotpos[*,0], thick=2
       oplot, wavesol_log[okHalog], cont_fit[okHalog], col=fsc_color('red'), thick=3
       oplot, wavesol_log[okHalog], best_fit[okHalog], col=fsc_color('blue'), thick=3
       oplot, wavesol_log[okHalog], line_fit[okHalog], col=fsc_color('green'), thick=3
       oplot, wavesol_lin[okHalin], (spec-cont_fit_lin)[okHalin]+10, col=fsc_color('magenta'), thick=2
       
       plot, wavesol_lin[okHblin], spec[okHblin], col=fsc_color('white'), position=plotpos[*,1], thick=2, /noerase
       oplot, wavesol_log[okHblog], cont_fit[okHblog], col=fsc_color('red'), thick=3
       oplot, wavesol_log[okHblog], best_fit[okHblog], col=fsc_color('blue'), thick=3
       oplot, wavesol_log[okHblog], line_fit[okHblog], col=fsc_color('green'), thick=3
       oplot, wavesol_lin[okHblin], (spec-cont_fit_lin)[okHblin]+10, col=fsc_color('magenta'), thick=2
       
       wait, 0.0001
       datacube[col,row,*] -= cont_fit_lin
    endif
  endfor
endfor

print, 'Saving the results...'
mwrfits, 0, outfile, /create
mwrfits, float(datacube), outfile, datahead
mwrfits, float(statcube), outfile, errhead

end
