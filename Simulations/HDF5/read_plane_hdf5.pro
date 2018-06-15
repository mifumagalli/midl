function read_plane_hdf5,infile,igz=igz,double=double
;
file_id    = h5f_open(infile)

; read the header data
group_id=h5g_open(file_id,'HEADER/')
att_id  = h5a_open_name(group_id,'nx')  & nx      = h5a_read(att_id) & h5a_close,att_id
att_id  = h5a_open_name(group_id,'ny')  & ny      = h5a_read(att_id) & h5a_close,att_id
att_id  = h5a_open_name(group_id,'nz')  & nz      = h5a_read(att_id) & h5a_close,att_id
att_id  = h5a_open_name(group_id,'iz')  & gz      = h5a_read(att_id) & h5a_close,att_id
h5g_close,group_id
; check whether this is right plane
if keyword_set(igz) then begin
    if (igz ne gz+1) then begin
        print,' wrong plane@ igz = ',igz,' gz = ',gz
    endif
endif
;
if keyword_set(double)then begin
    plane = dblarr(nx,ny)
endif else begin
    plane = fltarr(nx,ny)
endelse
;
grp_id   = h5g_open(file_id,'Data')
var_name =string('global_plane')
data_id = h5d_open(grp_id,var_name)
tmp     = h5d_read(data_id)
h5d_close,data_id
h5g_close,grp_id
plane   = tmp
;
h5f_close,file_id
;
return,plane
end
