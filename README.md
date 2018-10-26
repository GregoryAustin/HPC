USED gcc 7 and newest version of mpicc

To run the mpi 

mpirun -n cores mpi.o -i example1.txt -o outputmpi.txt

To run the omp 

./omp.o -i example1.txt -o outputomp.txt

To run the serial 

./serial.o -i example1.txt -o outputserial.txt 

NOTE: 
On my system the read from file works fine for reading the dcd file, but for the hex I needed to use a different read in system. 

This works on hex: 
int c = 0; 

dcdfile = malloc(sizeof(char*) * strlen(buf));

while (buf[c] != '\n') {
	dcdfile[c] = buf[c];
	c++;
}

dcdfile[c] = '\0';

This works on my machine: 
int c = 0; 

dcdfile = malloc(sizeof(char*) * strlen(buf));
while (buf[c] != '\0' && buf[c] != '\n') {
	dcdfile[c] = buf[c];
	c++;
}
// this gets rid of a line break before the \0 character
dcdfile[c-1] = '\0';
