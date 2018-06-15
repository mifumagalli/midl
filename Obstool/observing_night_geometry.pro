;+
;PURPOSE
;	to compute astrometric shifts from object1 to object2
;SYNTAX
;	observing_night_geometry, ra1, dec1, ra2, dec2
;INPUTS
;	ra1: right ascension in decimal degrees for object 1
;	dec1: declination in decimal degrees for object 1
;	ra2: right ascension in decimal degrees for object 2
;	dec2: declination in decimal degrees for object 2
;KEYWORDS
;	/sexi: set if your inputs are 3 element array 
;              in sexigesimal coordinates
;       /silent turn off the output
;OUTPUTS
;        
;		 
;Written by R. da Silva, UCSC
;MF add a couple of things (silent and output)
;preventing sleep deprived mistakes since 1-28-10
;-

PRO observing_night_geometry, ra1, dec1, ra2, dec2, sexi=sexi, silent=silent, $
                              pa=pa, shift=shift


if keyword_set(sexi) then begin
ra1=ten(ra1)*15D
dec1=ten(dec1)
ra2=ten(ra2)*15D
dec2=ten(dec2)
endif

ra_shift=(ra1-ra2)*3600.*cos(dec1*!dtor)
dec_shift=(dec1-dec2)*3600
if ra_shift GT 0 then ra_out2=' east' else ra_out2=' west'
if dec_shift GT 0 then dec_out2=' north' else dec_out2=' south'
ra_out=string(abs(ra_shift))+' arcseconds'+ra_out2
dec_out=string(abs(dec_shift))+' arcseconds'+dec_out2
posang, 1, ra1, dec1, ra2, dec2, pa
distance=djs_diff_angle(ra1, dec1, ra2, dec2)

posang, 1, ra1/15., dec1, ra2/15., dec2, pa
distance=djs_diff_angle(ra1, dec1, ra2, dec2)

if ~keyword_set(silent) then begin

    splog, 'PA = ', string(pa, format='(f8.3)'), ' , ',$
      string(pa+180, format='(f8.3)'),' , ', string(pa+360, format='(f8.3)')
    
    splog, 'These objects are ',distance*3600., ' arcseconds apart'
    
    splog, '-----------------------------------'
    splog, ' Shifts to move from obj1 to obj2:'
    splog, '-----------------------------------'
    splog, ra_out
    splog, dec_out
    splog, '-----------------------------------'
    splog, '    Move the telescope in this '
    splog, '  direction. I thought about this'
    splog, '           while awake'
    splog, '-----------------------------------'

endif


shift=[ra_shift,dec_shift]


end
