;+
;
; This is a procedure that takes an object structure and 
; produces a catalogue that can be used as input in the    
; bpz photoz code.
;
;
;  objstr  -->   name of the object structure with photometry
;  path    -->   where the object structure is
;
;
;-



pro bpz_gencat, objstr, path=path

;;open catalogue
obj=mrdfits(path+objstr,1)
nfilt=n_elements(obj[0].number)
nobj=n_elements(obj.number[0])


;;open out catalogue
bpz='bpz_'+strmid(objstr,0,strpos(objstr,'.fits'))+'.dat'

;;find zspec (distinguish between 0 and z)
obj.sdss_zspec=mk_finite(obj.sdss_zspec)
zspec=max(obj.sdss_zspec,dimension=1)
line=strarr(nobj)

openw, lun, bpz, /get_lun

for f=0, nfilt-1 do begin

    ;;calibrate color mag
    magcal=mk_finite(-2.5*alog10(obj.color_flux[f])+$
                     obj.zp_mag[f]-obj.color_term[f]-obj.am_term[f]-obj.galact_a[f])
    ;;compute error working backward rms calibration 
    cal_error_square=obj.errtot_magcal_cor[f]^2-(2.5*alog10(1.+obj.errtot_flux_cor[f]/obj.tot_flux_cor[f]))^2
    errcal=mk_finite(sqrt((2.5*alog10(1.+obj.errcolor_flux[f]/obj.color_flux[f]))^2+cal_error_square))
    ;;Set to 99 non detection (bpz makes distinction between non
    ;;detection 99 and non observed -99)
    uplim=where(magcal le 0 or errcal le 0 or magcal gt  40 or errcal gt 40, nul)
    magcal[uplim]=99
    errcal[uplim]=99

    splog, 'Calibrated magnitude for ', obj[0].filter[f]
    
    for o=0L, nobj-1 do begin
        
        ;;generate header
        if(o eq 0 and f eq 0) then  head='#id   z_spec   F'+obj[o].filter[f]+' E'+obj[o].filter[f]
        if(o eq 0 and f gt 0) then  head=head+' F'+obj[o].filter[f]+' E'+obj[o].filter[f]
        if(o eq 0 and f eq nfilt-1) then printf, lun, head
        
        ;;print stuff
        if(f eq 0) then line[o]=string(obj[o].number[0],zspec[o],magcal[o],errcal[o],format='(I8," ",F," ",F," ",F," ")')
        if(f gt 0) then line[o]=line[o]+string(magcal[o],errcal[o],format='(" ",F," ",F," ")') 
        if(f eq nfilt-1) then printf, lun, line[o]   
    
    endfor
    
    splog, 'Done filter ', obj[0].filter[f]
    
endfor

close, /all

end
