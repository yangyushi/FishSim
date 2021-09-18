prefix = /usr/local
CC = clang++
CPPFLAGS = -Wall -O3 -std=c++11 -I.. -I. -I$(prefix)/include -fPIC\
		   -Xpreprocessor -fopenmp -I/usr/local/opt/libomp/include -I$(prefix)/include/eigen3 
LINKFLAGS = -lomp


couzin: test_couzin.cpp couzin.o neighbour_list.o core.o
	$(CC) $(CPPFLAGS) $(LINKFLAGS) test_couzin.cpp couzin.o neighbour_list.o core.o
	rm -f test_couzin.xyz
	./a.out
	rm a.out

vicsek: test_vicsek.cpp vicsek.o neighbour_list.o network.o core.o
	$(CC) $(CPPFLAGS) $(LINKFLAGS) test_vicsek.cpp vicsek.o neighbour_list.o network.o core.o
	rm -f test.xyz
	./a.out
	rm a.out

network: test_network.cpp vicsek.o neighbour_list.o network.o core.o
	$(CC) $(CPPFLAGS) $(LINKFLAGS) test_network.cpp vicsek.o neighbour_list.o network.o core.o
	./a.out
	rm a.out

couzin.o: couzin.cpp
	$(CC) $(CPPFLAGS) -c couzin.cpp -o couzin.o

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
