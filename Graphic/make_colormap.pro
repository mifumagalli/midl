;+
;procedure che fa il ps di una simulation e aggiunge color bar
;range--> range of colorbar
;zoom is the the factor to zoom in a region 
;eps for encapsulated
;ps is the plate scale
;print save output in ps
;log for log scale
;zscale --> set zscale, if not set defualt if min max
;text allows label
;tlocation --> where to put the label
;colormap  allows to change the colormap
;grey     -->set to colormap 0
;reverse  reverse color map
;format  set the format for the colorbar ('F4.3')
;title -- display title
;coltitle --> title for the color bar
;colpos --> position for the color bar
;rotate --> 0-7 angle parameter for  rotation
;smooth --> the width of the smmothing function
;fits--> if set, im1 is passed as an array
;nocolbar  --> disable colorbar 
;x/ytitle   --> put a title in the bottom/left of the plot
;white     --> if set, the text is printed in black and background is
;              in white
;left      --> move the label to the left
;black     --> force black background
;circ      --> circle or radius pix, centered
;
;-

PRO make_colormap, img1, nameplot, RANGE=range, ZOOM=zoom, EPS=eps,$
                   PS=ps, print=print, log=log, zscale=zscale,  units=units, text=text, tlocation=tlocation,$
                   colormap=colormap,grey=grey,  reverse=reverse, format=format, title=title,$
                   coltitle=coltitle, rotate=rotate,$
                   smooth=smooth, fits=fits, colpos=colpos, nocolbar=nocolbar, xtitle=xtitle, ytitle=ytitle,$
                   white=white, left=left, black=black, _extra=extra, annotatecolor=annotatecolor, circ=circ


  ;set color map
  if ~keyword_set(colormap) then colormap=33
  if keyword_set(grey) then colormap=0
  if ~keyword_set(format) then format='(F5.2)'
  if keyword_set(white) then colormap=39
  if ~keyword_set(annotatecolor) then annotatecolor='black'

  if ~keyword_set(colpos) then colpos=[0.96,0.05,1.0,0.92]
  if ~keyword_set(reverse) then ctload, colormap
  if keyword_set(reverse) then ctload, colormap, /reverse
  
  IF keyword_set(fits) then fits=img1 else fits=MRDFITS(img1)
  
  if keyword_set(rotate) then fits=rotate(fits,rotate)
  if keyword_set(smooth) then fits=smooth(fits,smooth)

  
  ;;set zoom
  IF KEYWORD_SET(ZOOM) THEN BEGIN
      fsz=size(fits)
      xsize=fsz[1]-1
      ysize=fsz[2]-1
      DeltaX=(xsize-xsize/zoom)*0.5
      DeltaY=(ysize-ysize/zoom)*0.5
      fits=fits[DeltaX:(xsize-DeltaX),DeltaY:(ysize-DeltaY)]
      splog, 'Final size ', (size(fits))[1], (size(fits))[2]
      splog, 'Final size ', ps*(size(fits))[1], ps*(size(fits))[2]
      
  ENDIF
    
  IF keyword_set(print) THEN m_psopen, nameplot,/maxs, /encapsulated,/color
    
  ;;set scales
  IF KEYWORD_SET(range) THEN scale=range ELSE scale=[MIN(fits),MAX(fits)]
  IF KEYWORD_SET(zscale) THEN scale=ZSCALE_RANGE(fits)
 
  
  ;;if white set, saturate the background white
  IF KEYWORD_SET(white) THEN BEGIN
      ;;;fix white into red
      ;TVLCT, [[255], [0], [0]], 255
      ;;;put black into white
      ;TVLCT, [[255], [255], [255]], 0
      low=where(fits lt scale[0],nl)
      if(nl gt 0) then fits[low]=abs(scale[1]*10)
   ENDIF
  

  if KEYWORD_SET(black) THEN TVLCT, [[0], [0], [0]], 0
      

  !x.style=1
  !y.style=1
  
  IF KEYWORD_SET(log) THEN BEGIN
      IMDISP, ALOG10(fits), RANGE=ALOG10(scale), /axis, XTickformat='(A1)', YTickformat='(A1)',$
        XTICKLEN=1D-5,YTICKLEN=1D-5, title=title, xtitle=xtitle, ytitle=ytitle, $
        margin=0.025, charsize=2.5, position=[0.04,0.04,0.94,0.94], margin=0.
      if ~keyword_set(nocolbar) then colorbar, range=ALOG10(scale), FORMAT=format,$
        /vertical, title=coltitle, POSITION=colpos, annotatecolor=annotatecolor, _extra=extra
  ENDIF ELSE BEGIN
      IMDISP, fits, RANGE=scale, /axis, XTickformat='(A1)', YTickformat='(A1)',$
        XTICKLEN=1D-5,YTICKLEN=1D-5, title=title, xtitle=xtitle , ytitle=ytitle,$
        charsize=2.5, position=[0.01,0.02,0.9,0.92], margin=0.
      if ~keyword_set(nocolbar) then colorbar, range=scale,FORMAT=format, /vertical,$
        title=coltitle, POSITION=colpos, charsize=2., annotatecolor=annotatecolor, _extra=extra
  ENDELSE
  
 if keyword_set(left) then start=0.2 else start=0.8

 if keyword_set(white) then textcol=fsc_color("black") else textcol=fsc_color("white")


  IF KEYWORD_SET(units) THEN BEGIN
     xmax=N_ELEMENTS(fits[0,*])-1
     ymax=N_ELEMENTS(fits[0,*])-1
     oplot, [start*xmax-1./PS,start*xmax], [0.1*ymax,0.1*ymax], THICK=10, color=textcol
     lab='$Form.LABEL' 
     lab=units
     xyouts, (start*xmax-1./PS)+(0.5/PS), (0.05*ymax), lab, align=0.5, $
             CHARSIZE=1.5,CHARTHICK=3, color=textcol
  ENDIF
  
  
  if keyword_Set(circ) then x_oplotcirc, circ/ps, x0=((size(fits))[1])*0.5, y0=((size(fits))[2])*0.5, $
                                         color=textcol, line=2

  IF KEYWORD_SET(text) THEN BEGIN
     ;loadct, 0
     xyouts, tlocation[0],tlocation[1], text, CHARSIZE=1.5,CHARTHICK=2.5 , color=textcol 
  ENDIF
  
  IF keyword_set(print) THEN m_psclose




  ;;RESTORE COLOR BAR
  ctload, 39


END
