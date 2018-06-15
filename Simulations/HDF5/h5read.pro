function h5read,file,var
;
; modification history: written TT+RC May 10 2007
;
if (n_params() le 1) then begin
    print,' usage: h5read,file,var'
    print,' function to read variable "var" from hdf5 file "file"'
    return,0
endif
file_id    = h5f_open(file)
dataset_id = h5d_open(file_id,var)
data       = h5d_read(dataset_id)
h5d_close,dataset_id
h5f_close,file_id
return,data
end

