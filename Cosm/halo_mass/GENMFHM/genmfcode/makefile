#===============================================
# makefile for mass function program.
# Compile with real*8
#===============================================
genmf: 		genmf.o
	f77 -O genmf.o -o genmf
genmf.o:  genmf.f
	f77 -O -c genmf.f

