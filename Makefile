prefix=/usr/local

CPPFLAGS = -Wall -std=c++11 -Ofast
LINK = -I. -I${prefix}/include/eigen3 
LDLIBS = -L${prefix}/lib
EXTRAFLAGS = -fPIC -shared -mtune=native -march=native
PYFLAGS = `python3 -m pybind11 --includes` `python3-config --ldflags`

csimulate`python3-config --extension-suffix`: csimulate.cpp vicsek.o neighbour_list.o
	g++ $(CPPFLAGS) $(PYFLAGS) $(LINK) $(EXTRAFLAGS)\
		-o csimulate`python3-config --extension-suffix` csimulate.cpp vicsek.o neighbour_list.o
	rm -f test.xyz
	python3 simulate.py

vicsek.o: vicsek.cpp
	g++ $(CPPFLAGS) -c vicsek.cpp -o vicsek.o

neighbour_list.o: neighbour_list.cpp
	g++ $(CPPFLAGS) -c neighbour_list.cpp -o neighbour_list.o

test: main.cpp vicsek.o neighbour_list.o
	g++ $(LDLIBS) $(CPPFLAGS) $(LINK) main.cpp vicsek.o neighbour_list.o
	rm -f test.xyz
	./a.out
	rm a.out

clean:
	rm -f *.o
	rm -f a.out
	rm -f test.xyz
