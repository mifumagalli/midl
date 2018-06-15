;#############################################################################
; Generates an array of Gaussian emission lines to be used as templates in PPXF.
; logLam is the natural log of the wavelength of the templates in Angstrom.
; logLam should be the same as that of the stellar templates.
; FWHM_gal is the FWHM of the galaxy spectrum under study in Angstrom.
;
; V1.0: Michele Cappellari, Oxford, 7 January 2014
;#############################################################################
function ppxf_emission_lines, logLam_temp, FWHM_gal
;
; In this routine all lines are free to have independent intensities.
; One can fix the intensity ratio of different lines (e.g. the [OIII] doublet)
; by placing them in the same emission template
;
;        -----[OII]-----    Hdelta   Hgamma   Hbeta   -----[OIII]-----  [NI]    [OI]    -----[NII]-----   Halpha   -----[SII]-----
lines = [3726.03, 3728.82, 4101.76, 4340.47, 4861.33, 4958.92, 5006.84, 5199., 6300.30, 6548.03, 6583.41, 6562.80, 6716.47, 6730.85]
lam = exp(logLam_temp)
lines = lines[where((lines gt min(lam)) and (lines lt max(lam)))]
sigma = FWHM_gal/2.355 ; Assumes instrumental sigma is constant in Angstrom
emission_lines = dblarr(n_elements(logLam_temp),n_elements(lines))
for j=0,n_elements(lines)-1 do $
    emission_lines[*,j] = exp(-0.5d*((lam - lines[j])/sigma)^2)

return, emission_lines
end
;------------------------------------------------------------------------------
