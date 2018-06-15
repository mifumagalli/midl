function h5rd_multi,base,part_type,variable
;
; read multi-file snapshot variable
; input: base    : base name of file
; part_type=0..5  : particle type
; variable       : variable name
;

; read header information
inf           = base + '.0.hdf5'
NumPart_Total = h5ra(inf,'Header','NumPart_Total')
redshift      = h5ra(inf,'Header','Redshift')
nfiles        = h5ra(inf,'Header','NumFilesPerSnapshot')

; now read variable for each type
iread = 0l
for ifile=0,nfiles-1 do begin
   inf        = base + '.' + strtrim(string(ifile),1) + '.hdf5'   
   varname    = 'PartType' + strtrim(string(part_type),1) + '/' + strtrim(variable)
   data       = h5rd(inf,varname)
   NumPart_ThisFile = h5ra(inf,'Header','NumPart_ThisFile')
   print,' reading ',strtrim(inf),' np=',numpart_thisfile[part_type]
   if(ifile eq 0) then begin
      result = create_struct('z',redshift,'n',numpart_total[part_type],'var',replicate(data[0],NumPart_Total[part_type]))
   endif
   result.var[iread:iread+NumPart_ThisFile[part_type]-1] = data[*]
   iread = iread + NumPart_ThisFile[part_type]
endfor

return,result
end
