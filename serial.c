#define _GNU_SOURCE
#include <stdlib.h>
#include <stdio.h>
#include "includes/dcdplugin.c"
#include "priority_queue.c"
#include "linked_list.c"
#include <assert.h>
#include <time.h>
#include <omp.h>

unsigned int_size = sizeof(int); 

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
void createSet(int set_begin, int set_end, struct LNode** head_ref, int * count) {
    int i; 
    *count = 0; 
    for (i = set_end; i >= set_begin; --i) {
        // printf("PUSHING"); 
        push(head_ref, &i, int_size);
        *count += 1;  
    }
}

void calculateDistances3D(float *Ax, float *Ay, float *Az, int k, 
                          struct LNode* R, struct LNode* G,
                          Node **pq)
{
    int i;
    struct LNode * tempR = R; 
    struct LNode * tempG = G;
    // begin_a--; end_a--; begin_b--; end_b--; 
    while (tempR != NULL) {
        int rInt = *(int *) tempR->data;
        tempG = G; 
        while (tempG != NULL) {
            int bInt = *(int *) tempG->data;
            float dx = Ax[rInt] - Ax[bInt];
            float dy = Ay[rInt] - Ay[bInt];
            float dz = Az[rInt] - Az[bInt];

            float dx2 = dx * dx;
            float dy2 = dy * dy;
            float dz2 = dz * dz;

            float squaredDistance = dx2 + dy2 + dz2;

            // printf("%d %d %f\n", rInt, bInt, sqrtf(squaredDistance));

            pushQ(pq, rInt, bInt,  sqrtf(squaredDistance));

            tempG = tempG->next; 
        }

        tempR = tempR->next; 
    }
}



int main(int argc, char **argv) {
    /***************************************************************/
    /**********************Reading Input File **********************/
    /***************************************************************/
    double start_time, run_time;
    start_time = omp_get_wtime();

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
  	

  	struct LNode *setA = NULL; 
  	int countA = 0, countB = 0; 
    struct LNode *setB = NULL; 

    createSet(a_begin, a_end, &setA, &countA); // THESE ARE 100s 
    createSet(b_begin, b_end, &setB, &countB); // THESE ARE 100s 


    int natoms = 0;
    void *raw_data = open_dcd_read(dcdfile, "dcd", &natoms);//
    if (!raw_data) {
        printf("Please enter a valid name for the dcd file \n");//
        return 1;
    }
    dcdhandle *dcd = (dcdhandle *) raw_data;//

    molfile_timestep_t timestep;
    timestep.coords = (float *) malloc(3 * sizeof(float) * natoms);

    

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
        
        calculateDistances3D(dcd->x, dcd->y, dcd->z, 3, setA, setB, &pq);  


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
    
        printf("\n");


        // printf("Timestep %d\n", i);
        // printf("i: x    y      z\n");
        
        // if (i >= 10) break;
    }
    free(timestep.coords);    
    close_file_read(raw_data);

    run_time = omp_get_wtime() - start_time;
    printf("\n%lf seconds\n ",run_time);
    
    return 0;

}