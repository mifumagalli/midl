function h5wd,file,varname,var,create=create
                                ;
; modification history: written TT+RC May 10 2007
;
if (n_params() le 2) then begin
    print,' usage: h5w(file,varname,var)'
    print,' function to read variable "var" from hdf5 file "file"'
    return,0
endif

thisvar = size(var)
rank    = thisvar[0]
ndim    = intarr(rank)
if (rank gt 1) then begin
    for i=1,rank do ndim[i-1]=thisvar[i]
endif else begin
    ndim = thisvar[1]
endelse
vartype      = thisvar[rank+1]
dataspace_id = h5s_create_simple(ndim)

                                ; if file exists, open it, else create it
result   = file_info(file)
if (result.exists eq 1) then begin
    file_id     = h5f_open(file,/write)
    if keyword_set(create) then begin
        datatype_id  = H5T_IDL_CREATE(var)  
        dataset_id   = h5d_create(file_id,varname,datatype_id,dataspace_id)        
    endif else begin
        dataset_id  = h5d_open(file_id,varname)
    endelse
endif else begin
    file_id      = h5f_create(file)
    datatype_id  = H5T_IDL_CREATE(var)  
    dataset_id   = h5d_create(file_id,varname,datatype_id,dataspace_id)
endelse


; write the data
h5d_write,dataset_id,var

; close file and id
h5d_close,dataset_id
if(result.exists eq 1) then begin
    if keyword_set(create) then begin
        h5t_close,datatype_id        
    endif
endif else begin
    h5t_close,datatype_id
endelse
h5s_close,dataspace_id

h5f_close,file_id
;
return,0
end

