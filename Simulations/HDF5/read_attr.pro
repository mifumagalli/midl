function read_attr, file_name, group_name, attr_name, data_set=data_set

;
; reads attribute "attr_name" in group "group_name" of file "file_name"
; if keyword data_set is set, it is not a group but a data set
;

file_id = H5F_OPEN(file_name)
if keyword_set(data_set)then begin
    dataset_id = H5D_OPEN(file_id, group_name)
    attr_value = read_attribute_hdf5(dataset_id, attr_name)
    H5D_CLOSE, dataset_id
endif else begin
    group_id = H5G_OPEN(file_id, group_name)
    attr_value = read_attribute_hdf5(group_id, attr_name)
    H5G_CLOSE, group_id
endelse

H5F_CLOSE, file_id

return, attr_value

end
