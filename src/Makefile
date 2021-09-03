prefix = /usr/local
CC = clang++
CPPFLAGS = -Wall -O3 -std=c++11 -I.. -I. -I$(prefix)/include -fPIC\
		   -Xpreprocessor -fopenmp -I/usr/local/opt/libomp/include -I$(prefix)/include/eigen3 
LINKFLAGS = -lomp

main: main.cpp vicsek.o neighbour_list.o network.o
	$(CC) $(CPPFLAGS) $(LINKFLAGS) main.cpp vicsek.o neighbour_list.o network.o
	rm test.xyz
	./a.out

network.o: network.cpp
	$(CC) $(CPPFLAGS) -c network.cpp -o network.o

vicsek.o: vicsek.cpp
	$(CC) $(CPPFLAGS) -c vicsek.cpp -o vicsek.o

neighbour_list.o: neighbour_list.cpp
	$(CC) $(CPPFLAGS) -c neighbour_list.cpp -o neighbour_list.o

clean:
	rm -f *.o
	rm -f *.a
	rm -f a.out
	rm -f test.xyz
