;+
;
;Convert from angstrom to Hz
;
;
;
;
;
;-


pro ang2hz, ang, hz


hz=(2.9979D8)/(ang*1D-10)
print, hz

end
