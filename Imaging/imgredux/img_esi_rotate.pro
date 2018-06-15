;+
;
;
;procedure that for a given ROTATOR position, make the frame
;
;                  N
;                  ^
;                  |        
;                  |     
;            E <----
;
;
;
;-


pro img_esi_rotate, sci, wgt, rotpos, xlim, ylim, header=header

  ;;make this always positive
  pa_real=rotpos
  if(pa_real lt 0) then pa_real=360+pa_real
  
  ;;decide if image is rotated by multiple of 90 or 
  ;;by arbitrary angle 
  arbitrary=pa_real mod 90.

  if(arbitrary gt 2 and arbitrary lt 88) then begin
     
     ;;bring to closest integer of 90 and then repeat 
     if(arbitrary lt 45) then ang=-arbitrary else ang=90-arbitrary
     
     ;;15 deg is an empirical number used to tweak rotation
     sci=ROT(sci,ang-15,missing=0.)
     wgt=ROT(wgt,ang-15,missing=0.)

     pa_real=pa_real+ang
     
  endif
  
  ;;now do the 90 plus rotation
  
  ;;get the multiple of 90 you have to rotate for
  rot_ang=round(pa_real/90.)
  ;;make the rotation counterclockwise
  rot_ang=4-rot_ang
  ;;make the 360 angle a 0
  if(rot_ang eq 4) then rot_ang=0
  ;;take out an extra  90 degrees
  rot_ang=rot_ang+1
  if(rot_ang eq 4) then rot_ang=0
  ;;make the rotation
  sci=rotate(sci,rot_ang)
  wgt=rotate(wgt,rot_ang)
  ;;final reverse the X axis 
  sci=reverse(sci,2)  
  wgt=reverse(wgt,2)  
  
  ;;update the ccd limits if angle 1 and 3
  if(rot_ang eq 1 or rot_ang eq 3) then begin
     tmp=xlim
     xlim=ylim
     ylim=tmp
  endif
  
  ;;update header
  news=size(sci)
  sxaddpar, header, 'NAXIS1', news[1]
  sxaddpar, header, 'NAXIS2', news[2]
  
end

