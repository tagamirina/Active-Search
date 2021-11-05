CC       = g++
CFLAGS   = -O2
#CFLAGS   = -Wall

CVINC    = `pkg-config --cflags opencv`
CVLIB    = `pkg-config --libs opencv`
PATHS    = -I/usr/local/include -L/usr/local/lib -I/usr/local/include/opencv

clean :
	rm -f *.exe
	rm -f *.o
	rm -f *.a
	rm -f *.stackdump

all :
	make ${EXECS}

active : active.cpp
	${CC} ${CFLAGS} active.cpp ${CVINC} ${CVLIB} -o active.exe -g
	#g++ -O2 cptm.cpp `pkg-config --cflags opencv` `pkg-config --libs opencv` -o cptm.exe
