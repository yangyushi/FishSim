prefix = /usr/local
CC = g++
CPPFLAGS = -Wall -march=native -O3 -std=c++17 -fPIC \
		   -I. -I$(prefix)/include -I/usr/local/opt/libomp/include -I../extern/eigen3 \
		   -fopenmp

all: mcmc couzin_tank couzin vicsek network

test: mcmc couzin_tank couzin vicsek network
	./mcmc
	./vicsek
	./couzin
	./network
	./couzin_tank
	@rm *.xyz

mcmc: test_mc.cpp montecarlo.o core.o boundary.o
	$(CC) $(CPPFLAGS) test_mc.cpp  montecarlo.o core.o boundary.o -o mcmc

couzin_tank: test_couzin_tank.cpp couzin.o neighbour_list.o core.o boundary.o
	$(CC) $(CPPFLAGS)  test_couzin_tank.cpp couzin.o neighbour_list.o core.o boundary.o -o couzin_tank

couzin: test_couzin.cpp couzin.o neighbour_list.o core.o
	$(CC) $(CPPFLAGS) test_couzin.cpp couzin.o neighbour_list.o core.o boundary.o -o couzin

vicsek: test_vicsek.cpp vicsek.o neighbour_list.o network.o core.o
	$(CC) $(CPPFLAGS) test_vicsek.cpp vicsek.o neighbour_list.o network.o core.o -o vicsek

network: test_network.cpp vicsek.o neighbour_list.o network.o core.o
	$(CC) $(CPPFLAGS) test_network.cpp vicsek.o neighbour_list.o network.o core.o -o network

montecarlo.o: montecarlo.cpp
	$(CC) $(CPPFLAGS) -c montecarlo.cpp -o montecarlo.o

couzin.o: couzin.cpp boundary.o
	$(CC) $(CPPFLAGS) -c couzin.cpp -o couzin.o

boundary.o: boundary.cpp
	$(CC) $(CPPFLAGS) -c boundary.cpp -o boundary.o

network.o: network.cpp
	$(CC) $(CPPFLAGS) -c network.cpp -o network.o

vicsek.o: vicsek.cpp
	$(CC) $(CPPFLAGS) -c vicsek.cpp -o vicsek.o

neighbour_list.o: neighbour_list.cpp
	$(CC) $(CPPFLAGS) -c neighbour_list.cpp -o neighbour_list.o

core.o: core.cpp
	$(CC) $(CPPFLAGS) -c core.cpp -o core.o

clean:
	rm -f *.o
	rm -f *.a
	rm -f a.out
	rm -f test.xyz
	rm mcmc vicsek network couzin couzin_tank

.PHONY:
	clean test
