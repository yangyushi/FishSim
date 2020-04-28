#--------------------- Clang and openmp, working with my mac  ----------------
#CC = clang++
#CPPFLAGS = -Wall -O3 -std=c++11 -I.. -I. -I$(prefix)/include -fPIC\
		   -Xpreprocessor -fopenmp -lomp -I/usr/local/opt/libomp/include
#-----------------------------------------------------------------------------

#---------------------- GCC 9 and openmp, not working with my mac  -----------
CC = g++-9
CPPFLAGS = -Wall -O3 -std=c++11 -I.. -I. -I$(prefix)/include -fPIC -fopenmp
#-----------------------------------------------------------------------------

#--------------------- Clang , working with my mac  -------------------------
#CC = clang++
#CPPFLAGS = -Wall -O3 -std=c++11 -I.. -I. -I$(prefix)/include -fPIC
#-----------------------------------------------------------------------------

PYFLAGS = `python3 -m pybind11 --includes` `python3-config --ldflags` -L`python3-config --prefix`/lib 


csimulate`python3-config --extension-suffix`: csimulate.cpp libcsimulate.a
	$(CC) $(CPPFLAGS) $(PYFLAGS) -shared -o csimulate`python3-config --extension-suffix` csimulate.cpp libcsimulate.a
	mv *.so ../lib
	mv *.a ../lib

libcsimulate.a: vicsek.o neighbour_list.o
	ar rcs libcsimulate.a vicsek.o neighbour_list.o

vicsek.o: vicsek.cpp
	$(CC) $(CPPFLAGS) -c vicsek.cpp -o vicsek.o

neighbour_list.o: neighbour_list.cpp
	$(CC) $(CPPFLAGS) -c neighbour_list.cpp -o neighbour_list.o

sim: main.cpp vicsek.o neighbour_list.o
	$(CC) $(CPPFLAGS) main.cpp vicsek.o neighbour_list.o
	rm -f test.xyz
	time ./a.out
	rm a.out


test: test.cpp vicsek.o neighbour_list.o
	$(CC) $(CPPFLAGS) test.cpp vicsek.o neighbour_list.o
	time ./a.out
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
	rm -f *.a
	rm -f a.out
	rm -f test.xyz
