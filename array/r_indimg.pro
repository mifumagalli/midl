;+
;PURPOSE
;	to create a 2d image with values equal to its x index (or y index)
;SYNTAX
;	img=r_indimg(nx, ny, [ /y])
;INPUTS
;	nx: number of elements in x-direction
;	ny: number of elements in y-direciton
;KEYWORDS
;	/y: set if you want it to be the y index instead
;OPTIONAL SYNTAX
;	x=r_indimg(img, y=y)
;INPUTS
;	img: a 2d image that you want x and y index image for
;OUTPUTS
;	x: the x index image
;	y: the y index image
;~~~~~~~~~~~~~~~~~~ALTERNATIVE SYNTAX~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
;	ximg=r_indimg(img2, y=yimg)
;INPUTS		
;	img2: img that you want the index image of	
;OUTPUTS
;	ximg: x index image
;	yimg: y index image
;
;Written by R. da Silva, UCSC, 12-2-09
;Added optional syntax, R. da Silva, UCSC, 1-14-10
;-
FUNCTION r_indimg, nx, ny, y=y
if N_PARAMS() EQ 2 then begin
   if keyword_set(y) then return, findgen(ny)##replicate(1, nx) else $
	return, findgen(nx)#replicate(1, ny)
endif

if N_PARAMS() EQ 1 then begin
   dims=size(nx,/dim)
   y=findgen(dims[0])#replicate(1, dims[1])
   return, findgen(dims[1])##replicate(1, dims[0])
endif

end
