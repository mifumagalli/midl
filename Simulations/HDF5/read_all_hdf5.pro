; f128
base    = 'f128/save_g1_d2_c0_p'
nplanes = 128
f128    = fltarr(128,128,128)
for iz=0,nplanes-1 do begin $
  & file =base+string(iz,format='(i4.4)')+string('.d') $
  & f128[*,*,iz] = read_plane_hdf5(file) $
          & endfor

;
; f64
base    = 'f64/save_g1_d2_c0_p'
nplanes = 64
f64    = fltarr(64,64,64)
for iz=0,nplanes-1 do begin $
  & file =base+string(iz,format='(i4.4)')+string('.d') $
  & f64[*,*,iz] = read_plane_hdf5(file) $
          & endfor

; f32
base    = 'f32/save_g1_d2_c0_p'
nplanes = 32
f32    = fltarr(32,32,32)
for iz=0,nplanes-1 do begin $
  & file =base+string(iz,format='(i4.4)')+string('.d') $
  & f32[*,*,iz] = read_plane_hdf5(file) $
          & endfor

