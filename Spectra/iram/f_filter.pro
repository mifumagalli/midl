FUNCTION f_filter, input, A=A, B=B, tau=tau
;+
;PURPOSE
;	To fourier filter an offline spectrum to get the fitler shape
;SYNTAX
;	result=f_filter( input [tau=tau,A=A, B=B])
;INPUT
;	input: the spectra of off line measurement [1 d array]
;OUTPUT
;	filtered signal
;KEYWORDS
;	tau: how many points on each end to use
;	A: filtered inverse transform
;	B: unfiltered inverse transform
;
;NOTES
;	based on Carl Heiles's "How to Fourier Filter Your Reference Spectrum"
;	Need to optimize the use of tau. What value to use...
;-
nelem=n_elements(input)
if (n_elements(tau) EQ 0) then tau=32.
sym=[reverse(input[1:*]),input]
sym=[sym[0], sym]
A=fft(sym, /invers)
B=A
A[tau:n_elements(A)-tau-1]=0
output=fft(A)
output=output[nelem:*]
return, output
end
