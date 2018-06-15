;procedure che fa il ps di una singola immagine
;log for log scale
;lin for lin scale
;if log or lin non set, overrid with costum range 
;range is the scale range [min,max]
;zoom is the the factor to zoom in a region 
;ps is the plate scale
;if arcsec set, crea 
;text allows label
;eps for encapsulated

PRO image_simple, img1, nameplot, range, ZOOM=zoom, EPS=eps,$
    PS=ps, print=print, log=log, lin=lin, arcsec=arcsec, text=text


ctload, 0, /reverse


fits=MRDFITS(img1)


IF KEYWORD_SET(ZOOM) THEN BEGIN
xsize=N_ELEMENTS(fits[0,*])-1
ysize=N_ELEMENTS(fits[*,0])-1
DeltaX=(xsize-xsize/zoom)*0.5
DeltaY=(ysize-ysize/zoom)*0.5
fits=fits[DeltaX:(xsize-DeltaX),DeltaY:(ysize-DeltaY)]
ENDIF



IF keyword_set(print) THEN BEGIN
set_plot, "PS"
IF keyword_set(eps) THEN device, filename=nameplot, /encapsulated, /color, $
       ysize=10, xsize=10 ELSE  device, filename=nameplot, /color, $
       ysize=10, xsize=10
ENDIF


IF KEYWORD_SET(log) THEN scale=ZSCALE_RANGE(ALOG10(fits))
IF KEYWORD_SET(lin) THEN scale=ZSCALE_RANGE(fits)


IF KEYWORD_SET(log) THEN IMDISP, fits, RANGE=ALOG10(scale), /axis, XTickformat='(A1)', YTickformat='(A1)', XTICKLEN=1D-5,YTICKLEN=1D-5  
IF KEYWORD_SET(lin) THEN IMDISP, fits, RANGE=scale, /axis, XTickformat='(A1)', YTickformat='(A1)', XTICKLEN=1D-5,YTICKLEN=1D-5  
IMDISP, fits, RANGE=range, /axis, XTickformat='(A1)', YTickformat='(A1)', XTICKLEN=1D-5,YTICKLEN=1D-5  

IF KEYWORD_SET(arcsec) THEN BEGIN
xmax=N_ELEMENTS(fits[0,*])-1
ymax=N_ELEMENTS(fits[0,*])-1

loadct, 0

oplot, [0.8*xmax-2./PS,0.8*xmax], [0.1*ymax,0.1*ymax], THICK=2.5

lab='$Form.LABEL' 
lab='2"'
xyouts, (0.8*xmax-1.5/PS), (0.03*ymax), lab,CHARSIZE=1.5,CHARTHICK=2.5  
ENDIF


IF KEYWORD_SET(text) THEN BEGIN
;loadct, 39
xyouts, 0.55*xmax,0.5*ymax, "A", CHARSIZE=2.,CHARTHICK=2.5  
xyouts, 0.35*xmax,0.36*ymax, "B", CHARSIZE=2.,CHARTHICK=2.5  
ENDIF


IF keyword_set(print) THEN BEGIN
device, /close
set_plot, "X"
ENDIF

END
