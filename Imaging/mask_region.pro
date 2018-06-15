;+
;
;Read a ds9 region file with boxes and update a good pixel mask image.
;(0 to reject).
;
;regionfile    the ds9 region file, written in IRAF format.
;mask          input image to which the region mask is appended\no
;
;E.g. of one line of the region file:
;logical; box  2079.8144 4253.4935  4165.966 63.941079 3.9382058
;
;-


pro mask_region, regionfile, mask


;;read and parse the region file
readcol, regionfile, logical, box, xcenter, ycenter, xside, yside, ccrotation,$
  format='A,A,F,F,F,F,F'
nregion=n_elements(xcenter)


;;mask centers
xmas=fix(0.5*n_elements(mask[*,0]))
ymas=fix(0.5*n_elements(mask[0,*]))


;;loop over region and update mask
for reg=0, nregion-1 do begin
    
    ;;get edges
    x0=fix(xcenter[reg]-0.5*xside[reg])
    x1=fix(xcenter[reg]+0.5*xside[reg])
    y0=fix(ycenter[reg]-0.5*yside[reg])
    y1=fix(ycenter[reg]+0.5*yside[reg])
    

    ;;fix border
    if(x0 lt 0) then x0=0
    if(x1 ge n_elements(mask[*,0])) then x1=n_elements(mask[*,0])-1
    if(y0 lt 0) then y0=0
    if(y1 ge n_elements(mask[0,*])) then y1=n_elements(mask[0,*])-1
    
    
    ;;square box
    tmpmask=mask-mask+1
    tmpmask[x0:x1,y0:y1]=0.
    
    ;;rotate around center
    rotmask=ROT(tmpmask,360.-ccrotation[reg],1.,xcenter[reg],ycenter[reg],MISSING=1.)
    ;;now shift the centers
    shiftmask=shift(rotmask,xcenter[reg]-xmas,ycenter[reg]-ymas)
    ;;update mask
    mask=mask*shiftmask

    undefine, rotmask, tmpmask, shiftmask


endfor


end 
