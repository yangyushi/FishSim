CC = g++-9
CPPFLAGS = -std=c++14 -I. -Wall -fPIC
OPT = -O3
PYFLAGS = `python3 -m pybind11 --includes` `python3-config --ldflags` -L`python3-config --prefix`/lib -I/usr/local/include


csimulate`python3-config --extension-suffix`: csimulate.cpp vicsek.o neighbour_list.o
	$(CC) $(CPPFLAGS) $(OPT) $(PYFLAGS) -shared -o csimulate`python3-config --extension-suffix` csimulate.cpp vicsek.o neighbour_list.o
	rm -f test.xyz
	python3 simulate.py

vicsek.o: vicsek.cpp
	$(CC) $(CPPFLAGS) $(OPT) -c vicsek.cpp -o vicsek.o

neighbour_list.o: neighbour_list.cpp
	$(CC) $(CPPFLAGS) $(OPT) -c neighbour_list.cpp -o neighbour_list.o

test: main.cpp vicsek.o neighbour_list.o
	$(CC) $(CPPFLAGS) $(OPT) main.cpp vicsek.o neighbour_list.o
	rm -f test.xyz
	./a.out
	rm a.out

prof: main.cpp vicsek.o neighbour_list.o
	$(CC) $(CPPFLAGS) -O1 -pg main.cpp neighbour_list.cpp vicsek.cpp
	./a.out
	gprof ./a.out | gprof2dot | dot -Tpng -o output.png
	rm -f a.out
	rm -f test.xyz
	rm -f gmon.out

clean:
	rm -f *.o
	rm -f a.out
	rm -f test.xyz
