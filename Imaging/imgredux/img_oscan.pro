;Procedure that take the chip as input and subtract the overscan bias
;return a [xfull,yfull] data image bias subtracted 


pro img_oscan, chips, dataimage, layout=layout, instr=instr 
 
common ccdinfo, xfull, yfull, nchip, namp, biaslev, satur, goodata

;;get data sector min/max value
xmin=min(layout.data[0])
xmax=max(layout.data[1])
ymin=min(layout.y[0])
ymax=max(layout.y[1])

;;reconstruct the data image 
dataimage=make_array(xmax-xmin+1,ymax-ymin+1,/float)

;;run over each chip, debias each amplifier and trim the final image
;;to get rid of the bias section  

for amp=0, namp-1 do begin
    
   ;;data limit
   dxs=layout[amp].data[0]
   dxe=layout[amp].data[1]
   dys=layout[amp].y[0]
   dye=layout[amp].y[1]
   
   data=reform(chips[dxs:dxe,dys:dye])
   
   
   pxs=layout[amp].postdata[0]
   pxe=layout[amp].postdata[1]
   pys=layout[amp].y[0]
   pye=layout[amp].y[1]
   
   
   ;;exclude edges of overscan to avoid pattern in the bias
   if(instr eq 'FEINC') then overscan=reform(chips[pxs:pxe,pys:pye])
   if(instr ne 'FEINC' and instr ne 'LRIS') then overscan=reform(chips[pxs+5:pxe-5,pys:pye])
   if(instr eq 'LRIS') then begin
      ;;trim more aggressively to the right 
      overscan=reform(chips[pxs+10:pxe-5,pys:pye])
   endif
   
   ;;compute bias line by line
   linebias=djs_median(overscan,1)
   ;;smooth it
   smoo_pol=SAVGOL(30,30,0,3) 
   bias_smo=convol(linebias,smoo_pol,/EDGE_TRUNCATE)
   
   ;;make bias section
   bias=rebin(transpose(bias_smo),dxe-dxs+1,dye-dys+1)
   
   ;;compose trimmed mosaics
   xs = (dxs-xmin)
   xe = xs+(layout[amp].data[1]-layout[amp].data[0])
   ys = (dys-ymin)
   ye = ys+(layout[amp].y[1]-layout[amp].y[0])
   
   ;;subtract bias
   debias=data-bias
   dataimage[xs:xe,ys:ye]=debias
   
   ;;clean stuff
   undefine, debias, data, bias
   
endfor

if(instr eq 'LBC') then begin
   ;;reconstruct lbt mosaic
   mosaic=make_array(xfull,yfull)
   ;;fill in third chip
   mosaic[0:2047,0:4607]=dataimage[4096:6143,0:4607]
   ;;fill in second chip, leaving 18 pixx for chip gap
   mosaic[2122:4169,0:4607]=dataimage[2048:4095,0:4607]
   ;;fill in first chip, leaving 18 pixx for chip gap
   mosaic[4244:6291,0:4607]=dataimage[0:2047,0:4607]
   ;;rotate and fill in 4th chip,  leaving 18 pixy for chip gap
   chip4=rotate(dataimage[6144:8191,0:4607],1)
   mosaic[770:5377,4683:6730]=chip4
   
   ;;return
   dataimage=mosaic
   undefine, mosaic, chip4
endif

if(instr eq 'IMACS') then begin
   ;;reconstruct imacs mosaic
   mosaic=make_array(xfull,yfull)
     
     ;;;;;;;;;;;;;;;;
     ;; Chip order 
     ;; Note that top raw requires flip in XY 
     ;;  6 5 8 7  
     ;;  1 2 3 4
     ;; This puts a field in the south aligned with NE + 180 deg rot
     ;;;;;;;;;;;;;;;;
   
   order=[4,3,2,1,7,8,5,6]-1
   
   for amp=0, namp-1 do begin
    
      ;;data limit
      dxs=layout[order[amp]].data[0]
      dxe=layout[order[amp]].data[1]
      xs = (dxs-xmin)
      xe = xs+(layout[order[amp]].data[1]-layout[order[amp]].data[0])
      
      ;;set the y
      if(amp le 3) then begin
         
         ;;bottom 
         y0=ymin
         y1=ymax
         
         x0=8192-(amp+1)*2048
         x1=8191-amp*2048
         
         ;;bottom needs xy flip 
         mosaic[x0:x1,y0:y1]=dataimage[xs:xe,ymin:ymax]
         
      endif else begin

         ;;top 
         y0=ymin+ymax+1
         y1=ymax*2+1
         
         x0=8192-(amp-4+1)*2048
         x1=8191-(amp-4)*2048

         ;;top as it is
         mosaic[x0:x1,y0:y1]=rotate(dataimage[xs:xe,ymin:ymax],2)
         
      endelse
      
       
      
   endfor
   
   ;;return
   dataimage=mosaic
   undefine, mosaic

endif




                                ;plot, linebias, /ynozero
                                ;oplot, bias_smo, color=250
                                ;atv, dataimage, /block
                                ;stop
  
end
