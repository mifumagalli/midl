;+
;
;
;   This is a procedure that calls LABEL_REGION to find the
;   segmentation map of an image and perform basic stats of the regions
;
;   minval          the minimum level to consider 
;   maxval          the maximum level to consider
;   segm_map        in output, the segmentation map
;   all_neighbors   allows for 8-point neigh search
;   stats           a strucutre of basic info 
;
;-

pro find_connected, image, minval=minval, maxval=maxval, seg_map=seg_map, $
                    all_neighbors=all_neighbors, stats=stats

;;assign the default
if not keyword_set(maxval) then maxval=max(image)
if not keyword_set(minval) then minval=min(image)

;;create the mask image
imsiz=size(image)
mask=intarr(imsiz[1],imsiz[2])
inside=where(image ge minval and image le maxval,nfound)
if(nfound eq 0) then begin
    splog, 'No pixels inside range!'
    return
endif else mask[inside]=1

;;find cells
seg_map=LABEL_REGION(image*mask,all_neighbors=all_neighbors) 
;;free mask
undefine, mask

;;distance image
dist_ellipse, dist, [imsiz[1],imsiz[2]],floor(0.5*imsiz[1]),$
  floor(0.5*imsiz[2]),1.,0.

;;extract statistics of the image
nreg=max(seg_map)
npixel=fltarr(nreg)
meanrad=fltarr(nreg)
minrad=fltarr(nreg)
  
;;loop over regions
for i=1L, nreg do begin 
    ;;find pixel in current region
    inreg=where(seg_map eq i,npx)
    ;;set area
    npixel[i-1]=npx
    ;;find mean distance from center image
    meanrad[i-1]=total(dist[inreg]*image[inreg])/total(image[inreg])
    ;;find minimum distance from center
    minrad[i-1]=min(dist[inreg])
endfor
;;free dist
undefine, dist
stats={NUM_REG:1L*nreg,NUM_PIX:npixel,MEAN_DIST:meanrad,MIN_DIST:minrad}

end
