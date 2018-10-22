#define _GNU_SOURCE
#include <stdlib.h>
#include <stdio.h>
#include "../includes/dcdplugin.c"
#include "../priority_queue.c"
#include "../linked_list.c"
#include <assert.h>
#include <time.h>
#include <omp.h>

#define NUM_THREADS 4

unsigned int_size = sizeof(int); 

void print_array(int *array, int length)
{
    for (int i = 0; i < length; i++) { 
        printf("%d ", array[i]);

     }
     printf("\n");
}

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

/*
    Creates a LinkedList set given starting index and ending idnex  
*/
void createSet(int set_begin, int set_end, int list[]) {    

    int i; 
    for (i = 0; i <= (set_end - set_begin); i++) {
        list[i] = set_begin+i; 
    }
    
}

void calculateDistances3D(float *Ax, float *Ay, float *Az, int k, 
                          int listA[], int listB[], int aC, int bC, 
                          Node **pq)
{   
    // int nthreads  = 0;

    #pragma omp parallel
    {   
        Node * tQ = NULL; 
        // int nthrds = omp_get_num_threads();
        // if (omp_get_thread_num() == 0 ) nthreads = nthrds; 

        #pragma omp for collapse(2) 
        for (int i = 0; i < aC; ++i) {
            for (int j = 0; j < bC; ++j) {
                int aInd = listA[i]; 
                int bInd = listB[j];
                float dx = Ax[aInd] - Ax[bInd];
                float dy = Ay[aInd] - Ay[bInd];
                float dz = Az[aInd] - Az[bInd];

                float dx2 = dx * dx;
                float dy2 = dy * dy;
                float dz2 = dz * dz;

                float squaredDistance = dx2 + dy2 + dz2;
                squaredDistance = sqrtf(squaredDistance);

                pushQ(&tQ, aInd, bInd, squaredDistance, k);
            }
        }

        #pragma omp critical 
        while (!isEmpty(&tQ)) { 
            Node *pk = peek(&tQ);
            pushQ(pq, pk->a, pk->b, pk->priority, pk->k); 
            pop(&tQ); 
        } 
        
    }
    // printf("NUMBER OF THREADS: %d\n", nthreads); 
}

int main(int argc, char **argv) {

    double start_time, run_time;
    start_time = omp_get_wtime();
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
  	
    int aC = a_end - a_begin + 1;  
    int bC = b_end - b_begin + 1; 

  	int setA[aC]; 
    int setB[bC]; 
    
    createSet(a_begin, a_end, setA); 
    createSet(b_begin, b_end, setB); 

    omp_set_num_threads(NUM_THREADS); 

    int natoms = 0;
    void *raw_data = open_dcd_read(dcdfile, "dcd", &natoms);//
    if (!raw_data) {
        printf("Please enter a valid name for the dcd file \n");//
        return 1;
    }
    dcdhandle *dcd = (dcdhandle *) raw_data;//

    molfile_timestep_t timestep;
    timestep.coords = (float *) malloc(3 * sizeof(float) * natoms);


    /* TESTING TESTING TESTING */ 
    // Node *pq = NULL;

    // for (int i = 0; i < 20; i++) {
    //     pushQ(&pq, i, i, i, k-1);
    // }

    // while (!isEmpty(&pq)) { 
    //     Node *pk = peek(&pq);
    //     printf("%d, %d, %f\n", pk->a, pk->b, pk->priority); 
    //     pop(&pq); 
    //     count++;
    // } 


    /* TESTING TESTING TESTING */ 
    

    for (int i = 0; i < dcd->nsets; i++) {
        int rc = read_next_timestep(raw_data, natoms, &timestep);
        if (rc) {
            fprintf(stderr, "error in read_next_timestep on frame %d\n", i);
            return 1;
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
        int n = natoms; // > 10 ? 10 : natoms;
        
        Node *pq = NULL;
        

        // K-1 because that's the way the Queue works mates
        calculateDistances3D(dcd->x, dcd->y, dcd->z, k-1, setA, setB, aC, bC, &pq);  

        while (!isEmpty(&pq)) { 
            Node *pk = peek(&pq);
            // printf("%d, %d, %d, %f\n", i, pk->a, pk->b, pk->priority); 
            pop(&pq); 
        } 
    
        // printf("\n");


        // printf("Timestep %d\n", i);
        // printf("i: x    y      z\n");
        
        // if (i >= 8) break;
    }
    free(timestep.coords);    
    close_file_read(raw_data);

    run_time = omp_get_wtime() - start_time;
    printf("\n%lf seconds\n ",run_time);
    return 0;

}