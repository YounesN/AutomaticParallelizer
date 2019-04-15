CC=gcc-8
CFLAGS=-I.
DEPS = stencil.hpp
all: main.out generator.out

%.o: %.cpp $(DEPS)
	$(CC) -openmp -c -o $@ $< $(CFLAGS)

main.out: main.o stencil.o
	$(CC) -openmp -o main.out main.o stencil.o

generator.out: generator.o
	$(CC) -openmp -o generator.out generator.o

clean:
	rm -f *.out *.o *.out