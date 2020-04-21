CC = g++-9
CPPFLAGS = -O3 -std=c++14 -I. -Wall -fPIC
PYFLAGS = `python3 -m pybind11 --includes` `python3-config --ldflags` -L`python3-config --prefix`/lib -I/usr/local/include


csimulate`python3-config --extension-suffix`: csimulate.cpp vicsek.o neighbour_list.o
	$(CC) $(CPPFLAGS) $(PYFLAGS) -shared -o csimulate`python3-config --extension-suffix` csimulate.cpp vicsek.o neighbour_list.o
	rm -f test.xyz
	python3 simulate.py

vicsek.o: vicsek.cpp
	$(CC) $(CPPFLAGS) -c vicsek.cpp -o vicsek.o

neighbour_list.o: neighbour_list.cpp
	$(CC) $(CPPFLAGS) -c neighbour_list.cpp -o neighbour_list.o

test: main.cpp vicsek.o neighbour_list.o
	$(CC) $(CPPFLAGS) main.cpp vicsek.o neighbour_list.o
	rm -f test.xyz
	./a.out
	rm a.out

clean:
	rm -f *.o
	rm -f a.out
	rm -f test.xyz
