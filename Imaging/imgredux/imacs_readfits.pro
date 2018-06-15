;+
;;procedure that gets a imacs mosaic in inputs and makes a final fits
;;for science data. For test flats and focus frame, just disply the
;;image how it is and do not create any structure/mosaic.
;
;
;
;imagename  --> fits file as it comes from the telescope
;mosaic     --> if set to a variable, in output a single array full mosaic
;compose    --> generate final mosaic 
;
;written by MF Sept 2013
;-


function imacs_readfits, imagename, layout=layout, header=header, compose=compose
  
  
  ;;here some ccd specific as September 2009
  ;;relative position
  rel_dat=[0,2047]
  rel_pos=[2048,2111]
  rel_y=[0,4095]
  rel_fuy=4096
  rel_fux=16896
  
  ;;Derive absolute position - align all chips left to right 
  layout={DATA:[0,0],POSTDATA:[0,0],Y:[0,0]}
  layout=replicate(layout,8)
  for i=0L, 7 do begin
     layout[i].data=rel_dat+i*(rel_dat[1]-rel_dat[0]+1)
     layout[i].postdata=rel_pos+i*(rel_pos[1]-rel_pos[0]+1)+$
                        7*(rel_dat[1]-rel_dat[0]+1)
     layout[i].y=rel_y
     
  endfor
  
  ;;set mosaic
  mosaic=make_array(rel_fux,rel_fuy,/double)
  basename=strmid(imagename,0,strpos(imagename,"c1.fits.gz"))
  header=headfits(basename+"c1.fits.gz")
  
  for i=0, 7 do begin
     cc=mrdfits(basename+"c"+rstring(i+1)+".fits.gz",/fsc,/sil)
     mosaic[layout[i].data[0]:layout[i].data[1],0:rel_fuy-1]=cc[rel_dat[0]:rel_dat[1],0:4095]
     mosaic[layout[i].postdata[0]:layout[i].postdata[1],0:rel_fuy-1]=cc[rel_pos[0]:rel_pos[1],0:4095]
  endfor
  
  

  if keyword_set(compose) then begin
     
     ;;reconstruct imacs mosaic
     xfull=8192
     yfull=8192
     xmin=min(layout.data[0])
     xmax=max(layout.data[1])
     ymin=min(layout.y[0])
     ymax=max(layout.y[1])
     namp=8
     
     final=make_array(xfull,yfull)
     
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
           final[x0:x1,y0:y1]=mosaic[xs:xe,ymin:ymax]
           
        endif else begin
           
           ;;top 
           y0=ymin+ymax+1
           y1=ymax*2+1
           
           x0=8192-(amp-4+1)*2048
           x1=8191-(amp-4)*2048
           
           ;;top as it is
           final[x0:x1,y0:y1]=rotate(mosaic[xs:xe,ymin:ymax],2)
           
        endelse
        
     endfor

     xatv, final, /bl
     undefine, final
     
  endif
  

  return, mosaic
  
end
