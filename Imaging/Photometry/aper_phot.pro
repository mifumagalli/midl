;procedure to compute manual magnitude of an object in a FITS image.
;Consider this as a more advance XATV way to compute magnitudes,
;although it is not a fancy GUI.
;The image is assumed in e-/pix

;image       -->  the fits image 
;xpos        -->  x position of the object
;ypos        -->  y position of the object
;appix       -->  aperture in pixel 
;texp        -->  exposure time (used to measure the noise)
;Mags        -->  in output, the computed magnitude
;MagError    -->  in output, the error on the magnitude computed with
;                 standard signal-to-noise 


PRO aper_phot, image, xpos, ypos, appix, texp, Mags, MagError


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

if(x1 lt 0) then begin
xdiff=abs(x1)
endif
if(y1 lt 0) then begin
ydiff=abs(y1)
endif
if(x2 gt xlungo-1) then begin
xdiff=x2-xlungo+1
endif
if(y2 gt ylungo-1) then begin
ydiff=y2-ylungo+1
endif
 
diff=MAX([xdiff,ydiff]) 
 
boxsize=boxsize-diff

x1=xpos-boxsize
x2=xpos+boxsize
y1=ypos-boxsize 
y2=ypos+boxsize

 
cutfits=fits[x1:x2,y1:y2]


;empirical noise/mean ina box  
rad=40
pixels=xpix_circ(xpos,ypos,rad,/NOZERO, MAXX=xlungo, MAXY=ylungo, COUNT=PixMyaper)

;sigma clipping
sigma=50. 
fondo=0.
i=0


;assume gaussian noise (this can be modified to include correlated noise)
while (i lt 10) do begin
datatmp=double(fits[pixels[0,*],pixels[1,*]])
upper=fondo+3*sigma
lower=fondo-3*sigma
data=datatmp[where(datatmp gt lower and datatmp lt upper)]
skyhisto=histogram(data,min=min(data),binsize=0.002,locations=base)
yfit=gaussfit(base,skyhisto,parm,nterms=3)
sigma=parm[2]
fondo=parm[1] 
i=i+1
endwhile

;print, "Noise per pixel: ", sigma
;print, "Mean sky background: ", fondo


;sky subtract
fits=fits-fondo 


;find pixels in aperture and compute mag

Npix=dblarr(n_elements(appix))
counts=dblarr(n_elements(appix))

pixels=xpix_circ(xpos,ypos,appix,/nozero, maxx=xlungo, maxy=ylungo, count=number)
xindici=pixels[0,*]
yindici=pixels[1,*]
Npix=number
Counts=total(fits[xindici,yindici],/double)
mags=-2.5*alog10(counts)


;compute noise in aperture
noise=parm[2]*sqrt(npix)

SN=(Counts)/SQRT(noise^2+Counts/texp)

;compute error 
MagError=1.0857*SQRT(noise^2+Counts/texp)/Counts


;print, "Mag +/- Err : ", Mags, MagError 
;print, "Pixel: ", Npix
;print, "Signal to noise: ", SN



;some plot
 
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


END




