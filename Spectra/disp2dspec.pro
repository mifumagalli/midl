;this procedure displays a 2d spectrum that comes out of xidl
;long_reduce. 
;The wave image is reconstructed from one trace that oyu can specify.
;If not, it takes the brightest one.



;fitsimg  -> the reduce science frame
;trace    -> the trace to consider for the wave (number, not array index)
;imgout   -> in output, the image
;coadd2d  -> if set, deal with a 2d coadded image
;smooth   -> if set to a number, smooth by gaussian kernel

pro disp2dspec, fitsimg, trace, imgout=imgout, $
                coadd2d=coadd2d, smooth=smooth

if not keyword_set(coadd2d) then begin
;open image
    img=mrdfits(fitsimg,0)
    sky=mrdfits(fitsimg,2)

;open the structure
    str=mrdfits(fitsimg,5)
    
    if keyword_set(trace) then lambda=str[trace-1].WAVE_OPT $
    else begin
        max_flu=where(str.peakflux eq max(str.peakflux))
        lambda=str[max_flu].WAVE_OPT
    endelse
    
;make wave img
    wav=lambda ## replicate(1.,n_elements(img[*,0]))
    imgout=img-sky

endif else begin

    ;open image
    imgout=mrdfits(fitsimg,0)
    wav=mrdfits(fitsimg,4)

endelse

if keyword_set(smooth) then disp=smooth(imgout,smooth) else disp=imgout 

;displ
    
xatv, disp, wvimg=wav


end
