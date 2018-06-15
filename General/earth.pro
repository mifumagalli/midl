;+
;
;  plot the earth symbol
;
;-

pro earth 

A = FINDGEN(17)*(!PI*2/16.)
USERSYM, [COS(A), -1,1, 0,0, 0] , [SIN(A), 0, 0,0, -1, 1]

end
