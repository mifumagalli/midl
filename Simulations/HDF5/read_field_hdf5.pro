function read_field_hdf5,infile,double=double
;
file_id    = h5f_open(infile)

; read the header data
group_id=h5g_open(file_id,'HEADER/')
att_id  = h5a_open_name(group_id,'nx') & nx      = h5a_read(att_id) & h5a_close,att_id
att_id  = h5a_open_name(group_id,'ny') & ny      = h5a_read(att_id) & h5a_close,att_id
att_id  = h5a_open_name(group_id,'nz') & nz      = h5a_read(att_id) & h5a_close,att_id
h5g_close,group_id
;
data = dblarr(nx,ny,nz)
;
for iz=0,nz-1 do begin
    var_name =string('Planes/Plane')+string(iz,format='(i4.4)')
    data_id = h5d_open(file_id,var_name)
    plane   = h5d_read(data_id)
    h5d_close,data_id
    data[*,*,iz] = plane[*,*]
endfor
;
h5f_close,file_id
return,data
end
