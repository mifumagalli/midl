;+
; 
;
;  Procedure that parse SB99 output into a strucuture
;
;
;
;  For now, only quanta and spectrum are coded...
;
;
;
;
;
;-



pro sb99_parse, model, path=path

if ~keyword_set(path) then path='./'

;;make header model
sb99_mkhead, head, path+model+'.input1'

;;look for quanta

flag=FILE_TEST(path+model+'.quanta1')
if(flag eq 1) then sb99_quanta, quanta, path+model+'.quanta1', head else quanta=[0]


;;look for spectrum
flag=FILE_TEST(path+model+'.spectrum1')
if(flag eq '1') then sb99_spectrum, spectrum, path+model+'.spectrum1', head else spectrum=[0]


;;2-6 NOT CODED.........
;;>7not coded

;;write final file

;;Ext 0
mwrfits, [0], path+model+'.fits', head, /create
;;Ext 1 quanta
mwrfits, quanta, path+model+'.fits'
;;Ext 2
mwrfits, [0], path+ model+'.fits'
;;Ext 3
mwrfits, [0], path+model+'.fits'
;;Ext 4
mwrfits, [0], path+model+'.fits'
;;Ext 5
mwrfits, [0], path+model+'.fits'
;;Ext 6
mwrfits, [0], path+model+'.fits'
;;Ext 7
mwrfits, spectrum, path+model+'.fits'
;;Ext 8 ..... to be coded





end