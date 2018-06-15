;+
;PURPOSE
;	given a set of points regularly spaced in an x-array and a y-array
;	make a 2d image
;SYNTAX
;	img=reg_gridimg(zin, xin, yin [xout=xout, yout=yout, tol=tol])
;INPUT
;	zin: the value of a given point
;	xin: x coordinate of value
;	yin: y coordinate of value
;	tol: tolerance for matching values default to 1d-5
;OUPUTS
;	img: a 2d image spannig the regular values of xin and yin
;	xout: the x coordinate axis
;	yout the y coordinate axis
;NOTES:	
;	the input data MUST already be regularly gridded. This just
;	reorganizes the data
;Written by R. da Silva, UCSC, 11-30-09 
;-
FUNCTION reg_gridimg, zin, xin, yin, xout=xout, yout=yout, tol=tol
mmx=minmax(xin)
stepx=min(abs(xin[0]-xin[where(xin NE xin[0])]))
mmy=minmax(yin)
stepy=min(abs(yin[0]-yin[where(yin NE yin[0])]))
img=dblarr(round((mmx[1]-mmx[0])/stepx+1), round((mmy[1]-mmy[0])/stepy+1))

x1=0
for x=mmx[0], mmx[1], stepx do begin
   y1=0
   for y=mmy[0], mmy[1], stepy do begin
	img[x1, y1]=zin[where(neareq(xin, x, tol=tol) $
		AND neareq(yin, y, tol=tol))]
	y1++
   endfor
   x1++
endfor
xout=unique(xin)
yout=unique(yin)

return, img
end
