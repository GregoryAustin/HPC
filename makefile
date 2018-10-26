

all: $(EXES) 

serial$(EXE): serial-altalgo.c
	gcc -fopenmp serial-altalgo.c -o serial.o -lm

omp$(EXE): omp-altalgo.c
	gcc -fopenmp serial-altalgo.c -o omp.o -lm

mpi$(EXE): mpi_altalgo.c
	mpicc mpi_altalgo.c -o mpi.o -lm

clean: 
	rm *.o

.SUFFIXES:
.SUFFIXES: .c .cpp .$(OBJ)

.c.$(OBJ):
	$(CC) $(CFLAGS) -c $<

.cpp.$(OBJ):
	$(CC) $(CFLAGS) -c $<
