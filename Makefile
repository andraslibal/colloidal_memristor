XLIBS = -L/usr/X11R6/lib -lX11
XINCLUDE = -I/usr/X11R6/include/

memristor: src/memristor.c src/random_numrec.c src/util.c
	gcc src/memristor.c src/random_numrec.c src/util.c -o memristor -lm -O3 -mavx -ffast-math -freciprocal-math
directories:
	mkdir movies particles savedparams statistics
plot: plot.c 
	gcc plot.c -o plot -lm  $(XLIBS) $(XINCLUDE)
	
.PHONY: all swarm plot video parameter directories
all: swarm directories plot video parameter

clean:
	# rm -fr movies particles statistics savedparams
	rm -f memristor swarm.exe
