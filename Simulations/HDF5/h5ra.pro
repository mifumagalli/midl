function h5ra, file_name, group_name, attr_name, data_set=data_set
;
; modification history: written TT+RC May 10 2007
;
if (n_params() le 2) then begin
    print,' usage: h5ra(file_name,group_name,attribute_name, data_set=data_set_name)'
    print,' function to read attribute "attribute_name" from hdf5 file "file" in group "group_name"'
    print,'          if keyword data_set=data_set, attribuet is read for data_set'
    return,0
endif

file_id = H5F_OPEN(file_name)
if keyword_set(data_set)then begin
    dataset_id    = H5D_OPEN(file_id, group_name)
    attribute_id  = H5A_OPEN_NAME(dataset_id, attr_name)
    attr_value    = H5A_READ(attribute_id)
    H5A_CLOSE, attribute_id
    H5D_CLOSE, dataset_id
endif else begin
    group_id = H5G_OPEN(file_id, group_name)
    attribute_id  = H5A_OPEN_NAME(group_id, attr_name)
    attr_value    = H5A_READ(attribute_id)
    H5A_CLOSE, attribute_id
    H5G_CLOSE, group_id
endelse
H5F_CLOSE, file_id
return, attr_value
end
