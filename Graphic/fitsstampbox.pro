;+
;
; This procedure takes an image with wcs and make a postage 
; image around a given position, align NE. Default scaling is +/-5
; sigma bkg
;
;
; size      in arcsec
; range     manual range (min,max) of colorbar
; eps       for encapsulated
; print     save output in ps
; log       for log scale
; zscale    set zscale
; minmax    set minmax scale
; text      allows label
; tlocation where to put the label
; colormap  allows to change the colormap
; grey      set to colormap 0
; reverse   reverse color map
; format    set the format for the colorbar ('F4.3')
; title     display title
; coltitle  title for the color bar
; colpos    position for the color bar
; rotate    0-7 angle parameter for  rotation
; smooth    the width of the smoothing function
; nocolbar  disable colorbar 
; units     draw 1" line
; circ      put a circle of 10 pix radius around ra dec circ[0] circ[1] 
; nsig      the number of sky sigma to scale for
;-




pro fitsstampbox, img1, nameplot, ra, dec, size, range=range, eps=eps, $
                  print=print, log=log, zscale=zscale,  units=units, text=text, $
                  tlocation=tlocation,nsig=nsig,$
                  colormap=colormap, grey=grey,  reverse=reverse, format=format, $
                  title=title, coltitle=coltitle, $
                  rotate=rotate, smooth=smooth, colpos=colpos, nocolbar=nocolbar, $
                  circ=circ


;;set color map
if ~keyword_set(colormap) then colormap=39
if keyword_set(grey) then colormap=0
if ~keyword_set(reverse) then ctload, colormap
if keyword_set(reverse) then ctload, colormap, /reverse
if ~keyword_set(format) then format='(F5.2)'
if ~keyword_set(nsig) then nsig=5.


;;open image and header  
fits=mrdfits(img1,0,head,/sil)
   

;;align to NE
getrot, head, ne_rot, ne_cdelt
hrot,  fits, head, fits, head, ne_rot, -1, -1, 0


;;extract astrometric information
extast,head,astr
if keyword_set(rotate) then hrotate, fits, head, fits, head, rotate
if keyword_set(smooth) then fits=smooth(fits,smooth)


;;get the center position      
extast,head,astr
x_radec, ra, dec, rag, deg
ad2xy, rag, deg, astr, xcent, ycent
ps=astr.CD[1,1]*3600.
size_pix=size/ps

;;extract the box 
hextract, fits, head, fits, head, xcent-0.5*size_pix, xcent+0.5*size_pix,$
  ycent-0.5*size_pix, ycent+0.5*size_pix, /silent
      
if keyword_set(print) then m_psopen, nameplot, /land, /enc, /color
     


;;set default scales
mmm, fits, skymod, sky_sigma, sky_skew
fits=fits-skymod
scale=[-nsig*sky_sigma,nsig*sky_sigma]

;;overwrite is set
if keyword_set(range) then scale=range
if keyword_set(zscale) then scale=zscale_range(fits)
if keyword_set(minmax) then  scale=[min(fits),max(fits)]




!x.style=1
!y.style=1

  if keyword_set(log) then begin
      imdisp, alog10(fits), range=alog10(scale), /axis, xtickformat='(A1)',$
        ytickformat='(A1)', xticklen=1d-5,yticklen=1d-5, title=title 
      if ~keyword_set(nocolbar) then colorbar, range=alog10(scale),$
        format=format, /vertical, title=coltitle, position=colpos
  endif else begin
      imdisp, fits, range=scale, /axis, xtickformat='(A1)', $
        ytickformat='(A1)', xticklen=1d-5,yticklen=1d-5, title=title 
      if ~keyword_set(nocolbar) then colorbar, range=scale,$
        format=format, /vertical, title=coltitle, position=colpos
  endelse
  
  
  if keyword_set(units) then begin
     xmax=n_elements(fits[0,*])-1
     ymax=n_elements(fits[0,*])-1
     oplot, [0.8*xmax-1./ps,0.8*xmax], [0.1*ymax,0.1*ymax], thick=10, color=0
     lab=string(' 1" ')
     xyouts, (0.8*xmax-1./ps), (0.03*ymax), lab,charsize=1.5,charthick=2.5, color=0  
  endif
  
  
  if keyword_set(text) then begin
     loadct, 0
     xyouts, tlocation[0],tlocation[1], text, charsize=2.,charthick=2.5  
  endif
  
  ;;put circle
  if keyword_set(circ) then begin
      extast,head,astr
      x_radec, circ[0], circ[1], ragc, degc
      ad2xy, ragc, degc, astr, xc, yc
      tvcircle, 10., xc, yc, /data, color=0
  endif



  if keyword_set(print) then m_psclose
  
end
