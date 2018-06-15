function h5wa, file_name, group_name, attr_name, value, data_set=data_set
;
; modification history: written TT+RC May 10 2007
;
if (n_params() le 2) then begin
    print,' usage: h5wa(file_name,group_name,attribute_name, data_set=data_set_name)'
    print,' function to read attribute "attribute_name" from hdf5 file "file" in group "group_name"'
    print,'          if keyword data_set=data_set, attribuet is read for data_set'
endif

file_id = H5F_OPEN(file_name, /write)

; ; we begin by deleting the attribute
; if keyword_set(data_set)then begin
;     group_id      = H5D_OPEN(file_id, group_name)
;     attribute_id  = H5A_OPEN_NAME(group_id, attr_name)
;     h5a_delete, attribute_id, attr_name
;     H5D_CLOSE, group_id
; endif else begin
;     group_id = H5G_OPEN(file_id, group_name)
;     attribute_id  = H5A_OPEN_NAME(group_id, attr_name)
;     h5a_delete, group_id, attr_name
;     H5G_CLOSE, group_id
; endelse

; then we write the attribute
datatype_id  = H5T_IDL_CREATE(value)  
thisvar     = size(value)
rank        = thisvar[0]
ndim        = intarr(rank)
for i=1,rank do ndim[i-1]=thisvar[i]
dataspace_id = h5s_create_simple(ndim)

if keyword_set(data_set)then begin
    group_id      = H5D_OPEN(file_id, group_name)
    attribute_id  = H5A_CREATE(group_id, attr_name, datatype_id, dataspace_id)
    H5A_WRITE,attribute_id, value
    H5A_CLOSE, attribute_id
    H5D_CLOSE, group_id
endif else begin
    group_id = H5G_OPEN(file_id, group_name)
    attribute_id  = H5A_CREATE(group_id, attr_name, datatype_id, dataspace_id)
    H5A_WRITE,attribute_id, value
    H5A_CLOSE, attribute_id
    H5G_CLOSE, group_id
endelse
h5s_close,dataspace_id
h5t_close,datatype_id
H5F_CLOSE, file_id
;
return,0
end
