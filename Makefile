# Makefile for using PhilPlot
# Type: make blah.exe 

PHILSPLOT = /home/varnes/PhilsPlot

%.o : %.c
	gcc -c -I${PHILSPLOT}/include $< -o $@ -lm

%.exe : %.o 
	gcc -o $@ -I${PHILSPLOT}/include $< -L${PHILSPLOT}/lib -lphilsplot -lpgplot -lcpgplot -L/usr/X11R6/lib -lX11 -lpng -L/usr/lib/gcc/x86_64-redhat-linux/4.1.1 -lgfortran -lstdc++ -lm
