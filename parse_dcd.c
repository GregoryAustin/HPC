#define _GNU_SOURCE
#include <stdlib.h>
#include "includes/dcdplugin.c"

const char *dcd_file = "example_pn3_10RU_751frames.dcd";

/*
** Read and print first 10 timesteps and first 10 atoms per timestep in DCD files
** Look at dcdplugin.c and molfie_plugin.h for struct definitions
** The funtions called are very NOT thread safe, the following could cause data races:
**      - A file handle remains open from open_dcd_read till close_file_read
**      - dcd->x, dcd->y and dcd->z are overwritten each loop (among other properties)
**      - timestep.coords are overwritten each loop
**      - Probably other stuff? Also don't try to read the whole file into memory unless you have > 300GB of RAM
**
** If things start failing silently set NOISY to 1 on line 117 of dcdplugin.c (This might help or might not?)
*/
int main(int argc, char **argv) {
    int natoms = 0;
    void *raw_data = open_dcd_read(dcd_file, "dcd", &natoms);
    if (!raw_data) {
        printf("Please enter a valid name for the input file \n");
        return 1;
    }
    dcdhandle *dcd = (dcdhandle *) raw_data;

    molfile_timestep_t timestep;
    timestep.coords = (float *) malloc(3 * sizeof(float) * natoms);
    int read_failed = 0;
    for (int i = 0; i < dcd->nsets; ++i) {
        int rc = read_next_timestep(raw_data, natoms, &timestep);
        if (rc) {
            read_failed = 1;

            break;
        }
        /* At this point 
           dcd contains
           dcd->x = Array of X coordinates of all atoms for timestep i
           dcd->y = Array of Y coordinates of all atoms for timestep i
           dcd->z = Array of Z coordinates of all atoms for timestep i
           
           timestep contains
           timestep.coords = Array of packed XYZ coordinates of all atoms for timestep i
                             [X1, Y1, Z1, X2, Y2, Z2, ..., Xn, Yn, Zn] where n = natoms
          
           Both are overwritten next loop 
        */
        printf("Timestep %d\n", i);
        printf("i: x    y      z\n");
        int n = natoms > 10 ? 10 : natoms;
        for (int j = 0; j < n; ++j) {
            float *current = timestep.coords + 3 * j;
            printf("%d: %3.2f %3.2f %3.2f\n", j, current[0], current[1], current[2]);
        }
        printf("\n");
        if (i >= 10) break;
    }
    free(timestep.coords);    
    close_file_read(raw_data);
    if (read_failed) {
        return 2;
    }
    return 0;
}
