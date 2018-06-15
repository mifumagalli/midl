;+
;
;Take as an input x and y arrays and 
;generate an image of a 2d histo within 
;a given range of min/max
;  
;
; xarr       the array of x values
; yarr       the array of y value
; twod       the output image
; minx/maxx  min/max x-axis of output image
; miny/maxy  same for y axis
; binx/y     the pixel size for output image
; x/yaxis    in output, the array with the x/y values for the image edges.
;-




pro image_distribution,  xarr, yarr, twod, minx=minx, maxx=maxx, $
                         miny=miny, maxy=maxy, binx=binx, biny=biny, xaxis=xaxis, yaxis=yaxis
  

;;allocate image
nx=ceil((maxx-minx)/binx)
ny=ceil((maxy-miny)/biny)
twod=dblarr(nx,ny)


;;fill images with unit weight
value=replicate(1.,n_elements(xarr))

;;find indexes
xindx=(xarr-minx)/binx
yindx=(yarr-miny)/biny

;;populate
populate_image, twod, xindx, yindx, weight=value

;;output x/y values
xaxis=findgen(nx)*binx+minx
yaxis=findgen(ny)*biny+miny


end
