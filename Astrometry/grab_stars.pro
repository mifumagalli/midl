;grap the x and y position of a bunch of stars from a fits


pro grab_stars, fits, xref_star=xref_star,yref_star=yref_star, ps=ps, $
                dss=dss, scale=scale

  

   ;;-----------------------------------------------------
   ;;Provide a reference star position and do first shift
   ;;-----------------------------------------------------
    
  
  splog, 'Select your stars. Then, press right-button to continue...'

  
  ctload, 0, /reverse
  window, 0, xsize=1600, ysize=1600
  
  ;;reform and quick skysub
  imgzero=fits-djs_median(fits)
  if keyword_set(dss) then imdisp, bytscl(imgzero), /axis $
  else begin
     if keyword_set(scale) then $
        imdisp, imgzero, min=scale[0], max=scale[1], /axis $
     else begin
        sigma_sky=stddev(imgzero)
        imdisp, bytscl(imgzero, min=min(imgzero), max=min(imgzero)+sigma_sky), /axis
     endelse
  endelse
  
  ;;get first
  cursor, x, y, /wait, /data 
  ;;centroid
  cntrd, imgzero, x, y, xcen, ycen, 2./ps
  xref_tmp=xcen
  yref_tmp=ycen
  splog, 'Star Position ', xcen, ycen
  
  ;;display
  x_oplotcirc, 6./ps, x0=xcen, y0=ycen, /silent
  xyouts, xcen, ycen, "   1"
  
 
  ;;then continue
  for i=0, 100 do begin
     cursor, x, y, /down, /data 
     if(!mouse.button eq 1) then begin
        ;;centroid
        cntrd, imgzero, x, y, xcen, ycen, 2./ps
        xref_tmp=[xref_tmp,xcen]
        yref_tmp=[yref_tmp,ycen]
        splog, 'Star Position ', xcen, ycen
        
        ;;plot
        x_oplotcirc, 6./ps, x0=xcen, y0=ycen, /silent
        xyouts, xcen, ycen, string(i+2)
     endif
     if(!mouse.button eq 4) then break
  endfor
  
  undefine, imgzero
  
  ;;clean bad centroid
  goodcentr=where(xref_tmp GT 0 and yref_tmp GT 0,good)
  if(good gt 0) then begin
     xref_star=xref_tmp[goodcentr]
     yref_star=yref_tmp[goodcentr]
  endif
  
  ;;close windows
  wdelete, 0
  
end
