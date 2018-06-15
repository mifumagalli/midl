;+
;
;
; Build an image of an ellipse using the sextractor position parametrs
;
;
;
;-


pro sex_ellipse, distellipse, sizex=sizex, sizey=sizey, $
                 cxxwin_image=cxxwin_image, cyywin_image=cyywin_image, $ 
                 cxywin_image=cxywin_image, xwin_image=xwin_image, $
                 ywin_image=ywin_image

  ;;image x
  xarray=make_array(sizex,/index)
  ximage=rebin(xarray,sizex,sizey)
  ;;image y
  yarray=make_array(sizey,/index)
  yimage=rebin(transpose(yarray),sizex,sizey)
  
  ;;ellipse
  distellipse=CXXWIN_IMAGE*(ximage-XWIN_IMAGE)^2+$
              CYYWIN_IMAGE*(yimage-YWIN_IMAGE)^2+$
              CXYWIN_IMAGE*(ximage-XWIN_IMAGE)*$
              (yimage - YWIN_IMAGE)
  distellipse=sqrt(distellipse)/3.
  
end
