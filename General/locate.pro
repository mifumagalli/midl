;+
;
;
; Use the shell command locate to find an idl procedure that contains
; a name 
;
;
;
;
;-



pro locate, proc


spawn, 'locate '+proc+' | grep .pro', found
forprint, found, textout=1


end
