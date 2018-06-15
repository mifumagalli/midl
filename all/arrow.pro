; $Id: //depot/Release/ENVI50_IDL82/idl/idldir/lib/graphics/arrow.pro#2 $
;
; Copyright (c) 2005-2012, Exelis Visual Information Solutions, Inc. All
;       rights reserved. Unauthorized reproduction is prohibited.
;
;+
; :Description:
;    Create IDL Arrow graphic
;
; :Params:
;    parm1 : optional generic argument
;    parm2 : optional generic argument
;    parm3 : optional generic argument
;
; :Keywords:
;    ARROW_STYLE: 
;     0: ' --------'
;     1: ' ------->' --- default
;     2: ' <-------'
;     3: ' <------>'
;     4: ' >------>'
;     5: ' <------<'
;     
;     HEAD_INDENT:
;       Set this property to a floating-point value between -1 and +1 giving the
;       indentation of 
;       the back of the arrowhead along the shaft. A value of 0 gives a triangul
;       ar shape, 
;       a value of +1 will create an arrowhead that is just two lines, while a v
;        alue of -1 will 
;       create a diamond shape. The default is 0.4.
;       
;     HEAD_ANGLE:
;       Set this property to a floating-point value between 0 and 90 giving the 
;        angle in degrees 
;       of the arrowhead to the shaft. The default is 30.
;       
;    _REF_EXTRA
;
;-


;+
;
; Old idl procedure - function in >8.2?
;-

function arrow, arg1,y,z, _REF_EXTRA=_extra
  
  compile_opt idl2, hidden
  @graphic_error
  
  nparams = n_params()
  case (nparams) of
    1: begin
      return, call_function('arrow_internal_func', arg1, _EXTRA=_extra)
    end
    2: begin
      return, call_function('arrow_internal_func', arg1, y, _EXTRA=_extra)
    end
    3: begin
      return, call_function('arrow_internal_func', arg1, y, z, _EXTRA=_extra)
    end
    else: begin
      Message, 'Incorrect number of arguments.'
    end
  endcase 
end


;--------------------------------------------------------------------------
; This is the old ARROW procedure. We need to define its call here,
; and route the call to our internal .pro routine. Otherwise IDL will never
; find the old procedure name.
;
pro arrow, x0, y0, x1, y1, _REF_EXTRA=ex
  compile_opt hidden
  on_error, 2
  arrow_internal, x0, y0, x1, y1, _EXTRA=ex
end
