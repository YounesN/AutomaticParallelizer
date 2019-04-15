CC=g++
CFLAGS=-I.
DEPS = stencil.hpp
all: main.out generator.out

%.o: %.cpp $(DEPS)
	$(CC) -fopenmp -c -o $@ $< $(CFLAGS)

main.out: main.o stencil.o
	$(CC) -fopenmp -o main.out main.o stencil.o

generator.out: generator.o
	$(CC) -fopenmp -o generator.out generator.o

clean:
	rm -f *.out *.o *.out