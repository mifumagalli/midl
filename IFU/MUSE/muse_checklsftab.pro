;+
;
;  Check the line spread function out of muse pipe 
;
;-


pro muse_checklsftab, reduxpath=reduxpath


  if ~keyword_Set(reduxpath) then reduxpath=""

  ;;open output
  m_psopen, "lsf_check.ps", /land
  
  for i=0, 23 do begin

     ;;open table 
     tab=mrdfits(reduxpath+"LSF_TABLE-"+string(i+1,format='(I02)')+".fits",1,/sil)
     plothist, tab.lsf_width[0], bin=0.01, title='IFU '+string(i+1,format='(I02)'), $
               xtitle='GH Width'
     
     print, 'IFU ', i+1, ' ', median(tab.lsf_width[0])

  endfor
  
  m_psclose
  


end
