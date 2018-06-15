;+
;
; This procedure takes a fits file and resize it to smaller 
; dimension for presentation purpose.
; First it uses frebin to regrid. Then it makes an interpolation
; to make the image nicer.
;
;
;
; input  --> fits for input
; output --> fits resized to decent size for a ps
; nx,ny  --> the final size (defaul 300 x 300)
; scale  --> in output, two numbers for outsize/insize
; total  --> If set, values are added. This has to be used if the  
;            map is in flux units to preserve surface density.  
;            If not set, the values are averaged. This has to 
;            be used if the map is in units of surface density    
;            to preserve the final surface density.
; interp --> if set, after rebin do a cubic interpolation
;
;-



pro resizefits, input, output, nx=nx, ny=ny, scale=scale, total=total, $
                interp=interp
  
  ;;set default
  if ~keyword_set(nx) then nx=450
  if ~keyword_set(ny) then ny=450
  
  ;;grab the input size
  nxin=n_elements(input[*,0])
  nyin=n_elements(input[0,*])

  ;;compute scale factor 
  scale=[1.*nx/nxin,1.*ny/nyin]
  
  ;;allocate final grid 
  output=fltarr(nx,ny)
  
  ;;do interpolation on the new grid
  output=frebin(input,nx,ny,total=total) 

  
  if keyword_set(interp) then begin
     ;;apply the interpolation
     xax=make_array(nx,/index)
     yax=make_array(ny,/index)
     output=interpolate(output,xax,yax,cubic=-0.5,/grid)
  endif
  
  
end