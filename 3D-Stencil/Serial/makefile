CC=g++
CFLAGS=-I.
DEPS = stencil.hpp
all: main.out generator.out

%.o: %.cpp $(DEPS)
	$(CC) -c -o $@ $< $(CFLAGS)

main.out: main.o stencil.o
	$(CC) -o main.out main.o stencil.o

generator.out: generator.o
	$(CC) -o generator.out generator.o

clean:
	rm *.out *.o