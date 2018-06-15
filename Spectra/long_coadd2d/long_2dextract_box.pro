;Dummy version of long_2dextract that uses only one iteration
;and a simple box cart extraction.


;coadd_img  --> the coadded image from long_coadd2d
;xtrace     --> x position to be extracted
;radius     --> [optional] the width of the box
;extflux    --> output, the extracted flux
;extlambda  --> output, the extracted wavelenght 



pro long_2dextract_box, coadd_img, xtrace, radius=radius, $
                        extflux=extflux, extlambda=extlambda


;open long coadd stuff
sci=mrdfits(coadd_img,0)
lam=mrdfits(coadd_img,4)
mask=mrdfits(coadd_img,2)

;make ytrace 
ytrace=make_array(n_elements(sci[0,*]),/index)
xtrace=replicate(xtrace,n_elements(sci[0,*]))

;extract
extflux=extract_boxcar(sci*(mask EQ 0), xtrace, ytrace, radius=radius)
extlambda=extract_boxcar(lam, xtrace, ytrace, radius=radius)


splot, extlambda/(radius*2), extflux

xatv, mask


end
