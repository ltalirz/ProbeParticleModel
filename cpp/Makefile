CC=g++
CPPFLAGS=-std=c++11 -O3 -ffast-math -ftree-vectorize 
CPPLFLAGS= $(CPPFLAGS) -shared -Wl,-soname

all: PP GU MP STM

PP:
	$(CC) $(CPPFLAGS) -c -fPIC ProbeParticle.cpp 
	$(CC) $(CPPLFLAGS),ProbeParticle_lib.so -o ProbeParticle_lib.so ProbeParticle.o
GU:
	$(CC) $(CPPFLAGS) -c -fPIC GridUtils.cpp
	$(CC) $(CPPLFLAGS),GridUtils_lib.so -o GridUtils_lib.so GridUtils.o
MP:
	$(CC) $(CPPFLAGS) -c -fPIC Multipoles.cpp
	$(CC) $(CPPLFLAGS),Multipoles_lib.so -o Multipoles_lib.so Multipoles.o
STM:
	$(CC) $(CPPFLAGS) -c -fPIC ProbeSTM_0.0.10sp.cpp
	$(CC) $(CPPLFLAGS),ProbeSTM_0.0.10sp_lib.so -o ProbeSTM_0.0.10sp_lib.so ProbeSTM_0.0.10sp.o
clean:
	rm *.o *.so *.pyc
