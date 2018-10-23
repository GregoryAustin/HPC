#define _GNU_SOURCE
#include <stdlib.h>
#include <stdio.h>
#include "../includes/dcdplugin.c"
#include "../priority_queue.c"
#include "../linked_list.c"
#include <assert.h>
#include <time.h>
#include <omp.h>
#include "../dynam_arr.c"

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
void createSet(Array begins, Array ends, Array solos, struct LNode** head_ref) {
    // printf("Creating set... "); 

    for (int i = 0; i < begins.used; ++i) {
        for (int j = ends.array[i]; j >= begins.array[i]; --j) {
            push(head_ref, &j, int_size); 

        }
        // printf("Created a range"); 
    }

    for (int i = 0; i < solos.used; ++i) {
        push(head_ref, &(solos.array[i]), int_size); 
    }
}

// void createSet(Array begins, Array ends, Array solos, int list[]) {    

//     int count = 0; 
//     int size = 0; 

//     for (int i = 0; i < begins.used; ++i){ 
//         for (int j = 0; j <= (ends.array[i] - begins.array[i]); j++) {
//             list[size+j] = begins.array[i]+j; 
//         }
//         size += ends.array[i] - begins.array[i] + 1; 
//     }

//     for (int i = 0; i < solos.used; ++i) {
//         list[size+i] = solos.array[i];
//     }
    
// }


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

            pushQ(pq, rInt, bInt,  sqrtf(squaredDistance), k);

            tempG = tempG->next; 
        }

        tempR = tempR->next; 
    }
}


void getRangesAndItems(Array * begins, Array * ends, Array * solos, char ** tokens) {
    int cRange = 0;
    int items = 0;  
    int c = 0; 
    int begin, end;
    while (*(tokens+c) != NULL) {
        if (strpbrk(*(tokens+c), "-")) {
            
            char ** tokens2 = str_split(*(tokens+c), '-'); 

            sscanf(*(tokens2), "%d", &begin);
            sscanf(*(tokens2+1), "%d", &end);

            insertArray(begins, begin); 
            insertArray(ends, end); 
            cRange++; 
        } else {
            int temp; 
            // printf("Found no - in %s\n", *(tokens+c));
            sscanf(*(tokens+c), "%d", &temp);

            insertArray(solos, temp);  
            items++; 
        }

        c++; 
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
    int k;

    Array a_begins, b_begins; 
    Array a_ends, b_ends; 
    Array a_solos, b_solos; 

    initArray(&a_begins, 1); 
    initArray(&a_ends, 1); 
    initArray(&b_begins, 1); 
    initArray(&b_ends, 1); 
    initArray(&a_solos, 4); 
    initArray(&b_solos, 4); 


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
            tokens = str_split(buf, ',');
            getRangesAndItems(&a_begins, &a_ends, &a_solos, tokens); 
        } else if (count == 3) {
            tokens = str_split(buf, ',');
            getRangesAndItems(&b_begins, &b_ends, &b_solos, tokens); 
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


    printf("input: %s\nk: %d\n", dcdfile, k);


    printf("Set A: "); 

    for (int i = 0; i < a_begins.used; ++i) {
        printf("%d - %d, ", a_begins.array[i], a_ends.array[i]);
    }

    for (int i = 0; i < a_solos.used; ++i) {
        printf("%d, ", a_solos.array[i]); 
    }
    printf("\n"); 

    printf("Set B: "); 

    for (int i = 0; i < b_begins.used; ++i) {
        printf("%d - %d, ", b_begins.array[i], b_ends.array[i]);
    }

    for (int i = 0; i < b_solos.used; ++i) {
        printf("%d, ", b_solos.array[i]); 
    }
    printf("\n"); 

    fclose(ptr_file);
    /***************************************************************/
    /**************** Done Reading Input File **********************/
    /***************************************************************/
  	

  	struct LNode *setA = NULL; 
    struct LNode *setB = NULL; 

    // printf("about to create a set") ;

    createSet(a_begins, a_ends, a_solos, &setA); // THESE ARE 100s 
    createSet(b_begins, b_ends, b_solos, &setB); // THESE ARE 100s 

    // printList(setA, printInt); 

    // printList(setB, printInt); 


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
        
        calculateDistances3D(dcd->x, dcd->y, dcd->z, k-1, setA, setB, &pq);  

        // while (!isEmpty(&pq)) { 
        //     Node *pk = peek(&pq);
        //     printf("%d, %d, %d, %f\n", i, pk->a, pk->b, pk->priority); 
        //     pop(&pq); 
        // } 
    
        // printf("\n");


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