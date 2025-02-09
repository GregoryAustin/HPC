#define _GNU_SOURCE
#include <stdlib.h>
#include <stdio.h>
#include "includes/dcdplugin.c"
#include "includes/priority_queue.c"
#include "includes/linked_list.c"
#include <assert.h>
#include <time.h>
#include <omp.h>
#include "includes/mergeSort.c"
#include "includes/dynam_arr.c"
#include "includes/utility.c"

typedef int bool;
#define true 1
#define false 0

#define NUM_THREADS 4

unsigned int_size = sizeof(int); 


void getMinMax(float *axis, int * p1, int * p2, int fullSet[][2], int full) {
    float maxV = -9999999999.9;
    float minV = 9999999999.9;
    float temp;  

    for (int i = 0; i < full; ++i) {
        temp = axis[fullSet[i][0]]; 
        // printf("Temp: %f", temp);
        if (temp >= maxV) {
            maxV = temp; 
            *p2 = fullSet[i][0];
        } else if (temp <= minV) {
            minV = temp; 
            *p1 = fullSet[i][0];
        }
    }
} 

void calculateDistances3D(float *Ax, float *Ay, float *Az, int k, int fullSet[][2], int full, Node **pq) {
    float d1[full], d2[full], d3[full], d4[full], d5[full], d6[full], sum[full];
    int index[full];

    int p1, p2, p3, p4, p5, p6; 

    // biggest and smallest value sorted by x
    getMinMax(Ax, &p1, &p2,fullSet, full); 
    
    // biggest and smallest value sorted by y
    getMinMax(Ay, &p3, &p4, fullSet, full); 

    // biggest and smallest value sorted by z
    getMinMax(Az, &p5, &p6,fullSet, full); 

    //  Find distance of each point in p array from p1 and put its square in the d1 array. 

    #pragma omp parallel 
    {   
        #pragma omp for schedule(guided, 32) nowait
        for (int i = 0; i < full; ++i) 
            d1[i] = getSquaredDistance(Ax[p1], Ax[fullSet[i][0]], Ay[p1], Ay[fullSet[i][0]], Az[p1], Az[fullSet[i][0]]);

        //  Find distance of each point in p array from p2 and put its square in the d2 array. 
        #pragma omp for schedule(guided, 32) nowait
        for (int i = 0; i < full; ++i)
            d2[i] = getSquaredDistance(Ax[p2], Ax[fullSet[i][0]], Ay[p2], Ay[fullSet[i][0]], Az[p2], Az[fullSet[i][0]]);

        //  Find distance of each point in p array from p2 and put its square in the d2 array. 
        #pragma omp for schedule(guided, 32) nowait
        for (int i = 0; i < full; ++i) 
            d3[i] = getSquaredDistance(Ax[p3], Ax[fullSet[i][0]], Ay[p3], Ay[fullSet[i][0]], Az[p3], Az[fullSet[i][0]]);

        //  Find distance of each point in p array from p2 and put its square in the d2 array. 
        #pragma omp for schedule(guided, 32) nowait
        for (int i = 0; i < full; ++i) 
            d4[i] = getSquaredDistance(Ax[p4], Ax[fullSet[i][0]], Ay[p4], Ay[fullSet[i][0]], Az[p4], Az[fullSet[i][0]]);

        //  Find distance of each point in p array from p2 and put its square in the d2 array. 
        #pragma omp for schedule(guided, 32) nowait
        for (int i = 0; i < full; ++i) 
            d5[i] = getSquaredDistance(Ax[p5], Ax[fullSet[i][0]], Ay[p5], Ay[fullSet[i][0]], Az[p5], Az[fullSet[i][0]]);;
        
        //  Find distance of each point in p array from p2 and put its square in the d2 array. 
        #pragma omp for schedule(guided, 32) nowait
        for (int i = 0; i < full; ++i) 
            d6[i] = getSquaredDistance(Ax[p6], Ax[fullSet[i][0]], Ay[p6], Ay[fullSet[i][0]], Az[p6], Az[fullSet[i][0]]);

        // using a heuristic for the sum array 
        #pragma omp for schedule(guided, 32) nowait
        for (int i = 0; i < full; ++i){
            sum[i] = 11*d1[i] + 101*d2[i] + 547*d3[i] + 1009*d4[i] + 5501*d5[i] + 10007*d6[i];
            index[i] = i; 
        }
    }    
        // merge sorting sum array and swapping corresponding indexes 
        // #pragma omp task 
   mergeSort1(sum, index, 0, full-1); 
    
    int lookAhead = 100; 
    if (k != 0)
        lookAhead = 100 * sqrt(k); 

    // TODO: priority queue with fixed size 

    #pragma omp parallel 
    {
        Node * tQ = NULL; 

        #pragma omp for schedule(dynamic) nowait
        for (int i = 0; i < full; ++i) {
            for (int a = 1; a <= lookAhead && a+i < full; ++a) {
                int j = a+i; 

                int pointA = fullSet[index[i]][0];
                int pointB = fullSet[index[j]][0];  
                if(fullSet[index[i]][1] != fullSet[index[j]][1]) {
                    // printf("item p1: %d %d ; item p2: %d %d\n", fullSet[index[i]][0], fullSet[index[i]][1], fullSet[index[j]][0], fullSet[index[j]][1]);
                    float newPoint = getSquaredDistance(Ax[pointA], Ax[pointB], Ay[pointA], Ay[pointB], Az[pointA], Az[pointB]);

                    if (fullSet[index[i]][1] == 0) {
                        pushQ(&tQ, fullSet[index[i]][0], fullSet[index[j]][0], newPoint, k);
                    } else {
                        pushQ(&tQ, fullSet[index[j]][0], fullSet[index[i]][0], newPoint, k);
                    }
                }
            }
        }

        #pragma omp critical 
        while (!isEmpty(&tQ)) { 
            Node *pk = peek(&tQ);
            pushQ(pq, pk->a, pk->b, pk->priority, pk->k); 
            pop(&tQ); 
        } 
    }
    
    // printf("Smallest %d %d %lf", pA, pB, sqrtf(smallest));
    
}


int main(int argc, char **argv) {

    double start_time, run_time;
    start_time = omp_get_wtime();
    /***************************************************************/
    /**********************Reading Input File **********************/
    /***************************************************************/

    char *input_file; 
    char *output_file; 
    int nthreads = 0; 

    for(int i=0; i<argc; ++i)
    {   
        if (!strcmp(argv[i], "-i")) {
            if (i+1 < argc) 
                input_file = argv[i+1];
        } else if (!strcmp(argv[i], "-o")) {
            if (i+1 < argc) 
                output_file = argv[i+1];
        } else if (!strcmp(argv[i], "-n")) {
            if (i+1 < argc)
                nthreads = atoi(argv[i+1]);
        }
    }

    FILE *ptr_file; 
    char buf[1000];

    ptr_file = fopen(input_file, "r");

    FILE *output;
    output = fopen(output_file, "w");

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


    // printf("input: %s\nk: %d\n", dcdfile, k);

    // printf("Set A: "); 

    // for (int i = 0; i < a_begins.used; ++i) {
    //     printf("%d - %d, ", a_begins.array[i], a_ends.array[i]);
    // }

    // for (int i = 0; i < a_solos.used; ++i) {
    //     printf("%d, ", a_solos.array[i]); 
    // }
    // printf("\n"); 

    // printf("Set B: "); 

    // for (int i = 0; i < b_begins.used; ++i) {
    //     printf("%d - %d, ", b_begins.array[i], b_ends.array[i]);
    // }

    // for (int i = 0; i < b_solos.used; ++i) {
    //     printf("%d, ", b_solos.array[i]); 
    // }
    // printf("\n"); 

    fclose(ptr_file);
    /***************************************************************/
    /**************** Done Reading Input File **********************/
    /***************************************************************/
  	
    int aC = calculateSetSize(a_begins, a_ends, a_solos); 
    int bC = calculateSetSize(b_begins, b_ends, b_solos);

    // printf("aC size: %d\n", aC);
    // printf("bC size: %d\n", bC);

    int setA[aC][2]; 
    int setB[bC][2]; 

    
    // printf("About to create sets%d", fSize);;
    createSet(a_begins, a_ends, a_solos,0, setA); 
    createSet(b_begins, b_ends, b_solos,1, setB); 

    // print_array2(setB, bC); 
    
    int full[aC+bC][2];  //array_concat(setA, aC, setB, bC, int_size); 

    for (int i = 0; i < aC; i++) {
        full[i][0] = setA[i][0];
        full[i][1] = setA[i][1];
    }

    for (int i = 0; i < bC; i++) {
        full[aC+i][0] = setB[i][0];
        full[aC+i][1] = setB[i][1];
    }
    // printf("Combined sets");
    // print_array2(full, aC + bC);
    if (nthreads > 0)
        omp_set_num_threads(nthreads); 
    else 
        omp_set_num_threads(NUM_THREADS);

    int natoms = 0;
    void *raw_data = open_dcd_read(dcdfile, "dcd", &natoms);//
    if (!raw_data) {
        fprintf(stderr, "Please enter a valid name for the dcd file \n");//
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

        int n = natoms; // > 10 ? 10 : natoms;
        
        Node *pq = NULL; 
        calculateDistances3D(dcd->x, dcd->y, dcd->z, k-1, full, aC+bC, &pq);  

        while (!isEmpty(&pq)) { 
            Node *pk = peek(&pq);
            fprintf(output, "%d, %d, %d, %f\n", i, pk->a, pk->b, sqrtf(pk->priority)); 
            pop(&pq); 
        } 

        // cleaning up! 
        while (!isEmpty(&pq)) 
            pop(&pq); 
    }

    free(timestep.coords);    
    close_file_read(raw_data);

    run_time = omp_get_wtime() - start_time;
    fprintf(output, "\n%lf seconds\n ",run_time);
    fclose(output); 
    return 0;

}
