;procedure that for a given sdss spectrum return frequency and skysub
;flux



PRO get_sddsspec, sdssdata, lambda, sub_intens


spec=mrdfits(sdssdata,0,header)
cf0=SXPAR(header,"COEFF0")
cf1=SXPAR(header,"COEFF1")

index=fix(mkarr(0,N_ELEMENTS(spec[*,0]),1))

lambda=10^(cf0+cf1*index)
sub_intens=spec[*,0]

;units Flux   10-17 erg/cm2/s/A


END
