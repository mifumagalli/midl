;+
;
;
;procedure that for a given ROTATOR position, make the frame
;
;                  N
;                  ^
;                  |        This is tested for LRIS geometry 
;                  |        after the RED upgrade 
;            E <----        Write on 21/01/2010
;
;
;
;-


pro img_lris_rotate, sci, wgt, rotpos, xlim, ylim, header=header, noflip=noflip, rot_ang=rot_ang

;;take out the artificial -90
pa_real=rotpos+90
;;make this always positive
if(pa_real lt 0) then pa_real=360+pa_real
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

;;final reverse the X axis if new detector
date=fxpar(header,"DATE")

if(date gt '2009-04-25') and ~keyword_set(noflip) then begin  
   sci=reverse(sci)  
   wgt=reverse(wgt)  
endif

;;update the ccd limits if angle 1 and 3
if(rot_ang eq 1 or rot_ang eq 3) then begin
   tmp=xlim
   xlim=ylim
   ylim=tmp
endif


end

