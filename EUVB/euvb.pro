;+
;
;  Procedure that returns the EUVB spectrum at a given redshift 
;  (rounded to closest integer)
;
;
;
;
;-

pro euvb, wave, j22, red=red

 if not keyword_set(red) then red=3.0

 

 ;;load table 
 openr, lun, getenv('MIDL')+'/EUVB/HM2012_QG_uvb.dat', /get_lun

 ;;header
 head=strarr(20)
 readf, lun, head

 ;;redshift
 redshift=fltarr(60)
 readf, lun, redshift
 

 ;;wave and flux
 wave=fltarr(575)
 uvb=fltarr(575,60)

 for i=0, 574 do begin
   
    bufw=0.
    bufi=fltarr(60)
    readf, lun, bufw, bufi

    wave[i]=bufw
    uvb[i,*]=bufi

 endfor
 close, lun

 ;;find closest redshift
 zclose=min(redshift-red,zindx,/abs)+red
 j22=uvb[*,zindx]
 
end
