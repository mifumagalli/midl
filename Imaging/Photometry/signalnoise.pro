;this procedure computes the empirical SN in a region 
;over an image in e-/s.
;This is adding also the noise of the object which is in the aperture
;No correlation is taken into account. Use effectivenoise instead.
;
;image    --> the fits file
;xpos     --> x position 
;ypos     --> y position 
;appix    --> the aperture in pixel
;texp     --> the exposure time (used to compute the correct SN)




PRO signalnoise, image, xpos, ypos, appix, texp 

loadct, 0




;open fits file
fits=mrdfits(image,0,header,/fscale,/silent)
fits=fits 

;cut a box
xlungo=N_ELEMENTS(fits[*,0])
ylungo=N_ELEMENTS(fits[0,*])


boxsize=50.
x1=xpos-boxsize
x2=xpos+boxsize
y1=ypos-boxsize 
y2=ypos+boxsize

xdiff=0.
ydiff=0.

IF(x1 LT 0) THEN BEGIN
xdiff=ABS(x1)
ENDIF
IF(y1 LT 0) THEN BEGIN
ydiff=ABS(y1)
ENDIF
IF(x2 GT xlungo-1) THEN BEGIN
xdiff=x2-xlungo+1
ENDIF
IF(y2 GT ylungo-1) THEN BEGIN
ydiff=y2-ylungo+1
ENDIF
 
diff=MAX([xdiff,ydiff]) 
 
boxsize=boxsize-diff

x1=xpos-boxsize
x2=xpos+boxsize
y1=ypos-boxsize 
y2=ypos+boxsize

 
cutfits=fits[x1:x2,y1:y2]


;empirical noise/mean  
rad=40
pixels=xpix_circ(xpos,ypos,rad,/NOZERO, MAXX=xlungo, MAXY=ylungo, COUNT=PixMyaper)

;sigma clipping (here gaussianity is assumed, but it can be improved
;to account for non gaussian tails)
sigma=50. 
fondo=0.
i=0
WHILE (i LT 10) DO BEGIN
datatmp=double(fits[pixels[0,*],pixels[1,*]])
upper=fondo+3*sigma
lower=fondo-3*sigma
data=datatmp[where(datatmp GT lower AND datatmp LT upper)]
skyhisto=HISTOGRAM(data,MIN=MIN(data),BINSIZE=0.002,LOCATIONS=base)
YFIT=GAUSSFIT(base,skyhisto,PARM,NTERMS=3)
sigma=PARM[2]
fondo=PARM[1] 
i=i+1
ENDWHILE

print, "Noise per pixel: ", sigma
print, "Mean sky background: ", fondo


;sky sub
fits=fits-fondo 


;find pixels in aperture and compute mag

Npix=DBLARR(N_ELEMENTS(appix))
Counts=DBLARR(N_ELEMENTS(appix))

pixels=xpix_circ(xpos,ypos,appix,/NOZERO, MAXX=xlungo, MAXY=ylungo, COUNT=number)
xindici=pixels[0,*]
yindici=pixels[1,*]
Npix=number
Counts=TOTAL(fits[xindici,yindici],/double)


; compute noise in aperture
noise=PARM[2]*SQRT(Npix)

SN=(Counts)/SQRT(noise^2+Counts/texp)


print, "Counts: ", Counts
print, "Pixel: ", Npix
print, "Signal to noise: ", SN




;un po' di plot
 
!P.MULTI=[0,2,1]
Result = ZSCALE_RANGE(cutfits)
IMDISP, cutfits, RANGE=Result, /axis 

PLOTS, CIRCLE(boxsize, boxsize, appix), /DATA, thick=3, color=0
PLOTS, CIRCLE(boxsize, boxsize, rad), /DATA, thick=3, color=255

Z=double((base-PARM[1])/PARM[2])
GY=PARM[0]*EXP(-Z^2/2.)

plot, base, skyhisto, psym=10, xtitle="Sky Counts", ytitle="N" 
oplot, base, GY
oplot, [PARM[1], PARM[1]], [-1000,1000], line=2

!P.MULTI=0

usecolor

END




