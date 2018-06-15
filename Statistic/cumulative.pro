;+
;
;procedure that for a give array computes the cumulative distribution
;
; histo --> the input array
; cumulative --> the output cumulative distribution 
;                (sum of the various elements)
; bin--> the bin size of the differentail distribution
;
;-

PRO cumulative, histo, cumulative, BIN=bin
  
  if ~keyword_set(bin) then bin=1.
  
  dim_histo=n_elements(histo)
  cumulative=fltarr(dim_histo)
  for i=0, dim_histo-1 do cumulative[i]=total(histo[0:i]*bin)
  
end
