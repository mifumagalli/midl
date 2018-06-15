;+
;
;    This is a utility that figures out the location 
;    and number of lines grab a filter in 
;    master.FILTER.RES and master.FILTER.RES.info
;
;
;-



pro getfilter, filter_id, lambda=lambda, trans=trans

  
  ;;read master files info
  readcol, getenv('MIDL')+'/Photz/filter/master.FILTER.RES.info', nfilt, $
           startlines, lineinfilt, format='L,A,L', /sil
  
  ;;trim semicolum
  p=strpos(startlines,':')
  startline=nfilt-nfilt
  nold=n_elements(p)
  for i=0, nold-1 do startline[i]=1L*strmid(startlines[i],0,p[i])
  

  ;;find this values
  thistart=startline[filter_id-1]
  thisnline=lineinfilt[filter_id-1]
  
  
  ;;read the filter
  readcol, getenv('MIDL')+'/Photz/filter/master.FILTER.RES', id, lambda, trans, $
           skipline=thistart, numline=thisnline, /silent

  
end
