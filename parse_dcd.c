#define _GNU_SOURCE
#include <stdlib.h>
#include <stdio.h>
#include "includes/dcdplugin.c"
#include "priority_queue.c"
#include <assert.h>

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

char** str_split(char* a_str, const char a_delim)
{
    char** result    = 0;
    size_t count     = 0;
    char* tmp        = a_str;
    char* last_comma = 0;
    char delim[2];
    delim[0] = a_delim;
    delim[1] = 0;

    /* Count how many elements will be extracted. */
    while (*tmp)
    {
        if (a_delim == *tmp)
        {
            count++;
            last_comma = tmp;
        }
        tmp++;
    }

    /* Add space for trailing token. */
    count += last_comma < (a_str + strlen(a_str) - 1);

    /* Add space for terminating null string so caller
       knows where the list of returned strings ends. */
    count++;

    result = malloc(sizeof(char*) * count);

    if (result)
    {
        size_t idx  = 0;
        char* token = strtok(a_str, delim);

        while (token)
        {
            assert(idx < count);
            *(result + idx++) = strdup(token);
            token = strtok(0, delim);
        }
        assert(idx == count - 1);
        *(result + idx) = 0;
    }

    return result;
}

void calculateDistances3D(float *Ax, float *Ay, float *Az, int k, 
                          int begin_a, int end_a, int begin_b, int end_b,
                          Node **pq)
{
    int i;

    // begin_a--; end_a--; begin_b--; end_b--; 

    for (i = begin_a; i <= end_a; i++) {
        int j; 
        for (j = begin_b; j <= end_b; j++) {
            float dx = Ax[i] - Ax[j];
            float dy = Ay[i] - Ay[j];
            float dz = Az[i] - Az[j];

            float dx2 = dx * dx;
            float dy2 = dy * dy;
            float dz2 = dz * dz;

            float squaredDistance = dx2 + dy2 + dz2;

            push(pq, i, j, sqrtf(squaredDistance));


            // pairs[i] = (pair) {i, j, sqrtf(squaredDistance)};
        }
    }
}




int main(int argc, char **argv) {
    /***************************************************************/
    /**********************Reading Input File **********************/
    /***************************************************************/

    char *input_file; 
    char *output_file; 
    for(int i=0; i<argc; ++i)
    {   
        if (!strcmp(argv[i], "-i")) {
            if (i+1 < argc) 
                input_file = argv[i+1];
        } else if (!strcmp(argv[i], "-o")) {
            if (i+1 < argc) 
                output_file = argv[i+1];
        }
    }

    FILE *ptr_file; 
    char buf[1000];

    ptr_file = fopen(input_file, "r");

    if (!ptr_file)
        return 1;

    int count = 0; 
    int k, a_begin, a_end, b_begin, b_end;
    char **tokens; 
    char *dcdfile; 

    while (fgets(buf,1000, ptr_file)!=NULL) {
        if (count == 0) {
            int c = 0; 


            dcdfile = malloc(sizeof(char*) * strlen(buf));
             while (buf[c] != '\0' && buf[c] != '\n') {
                dcdfile[c] = buf[c];
                c++;
             }
             // this gets rid of a line break before the \0 character
             // TODO: a safer way to do this (ignore \n's) 
             dcdfile[c-1] = '\0';
        } else if (count == 1)
            sscanf(buf, "%d", &k);
        else if (count == 2) {
            //TODO: handle itemized list 
            //TODO: handle ranges within itemized list 
            tokens = str_split(buf, '-');
            sscanf (*(tokens), "%d", &a_begin);
            sscanf(*(tokens+1), "%d", &a_end);
        } else if (count == 3) {
            tokens = str_split(buf, '-');
            sscanf (*(tokens), "%d", &b_begin);
            sscanf(*(tokens+1), "%d", &b_end);
        }
        // a little bit of mem mgmt
        if (count == 2 || count == 3) {
            free(*(tokens));
            free(*(tokens+1));
            free(tokens);
        }
        count++; 
    }
    // printf("input: %s", dcdfile);


    printf("input: %s\nk: %d\nA: %d - %d\nB: %d - %d\n", dcdfile, k, a_begin, a_end, b_begin, b_end);

    fclose(ptr_file);
    /***************************************************************/
    /**************** Done Reading Input File **********************/
    /***************************************************************/
    
    int natoms = 0;
    void *raw_data = open_dcd_read(dcdfile, "dcd", &natoms);
    if (!raw_data) {
        printf("Please enter a valid name for the dcd file \n");
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


        // printf("Timestep %d\n", i);
        // printf("i: x    y      z\n");
        int n = natoms; // > 10 ? 10 : natoms;
        
        Node *pq = newNode(1,1,999999); 
        // printf("HERE");
        calculateDistances3D(dcd->x, dcd->y, dcd->z, 3, a_begin, a_end, b_begin, b_end, &pq);  


        int count = 0; 
        while (!isEmpty(&pq) && count < 3) { 
            Node *pk = peek(&pq);
            printf("%d, %d, %d, %f\n", i, pk->a, pk->b, pk->priority); 
            pop(&pq); 
            count++;
        } 

        // cleaning up! 
        while (!isEmpty(&pq)) 
            pop(&pq); 
        
      

        // for (int j = 0; j < n; ++j) {
        //     float *current = timestep.coords + 3 * j;
        //     printf("%d: %3.2f %3.2f %3.2f\n", j, current[0], current[1], current[2]);
        // }
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