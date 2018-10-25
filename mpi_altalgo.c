#define _GNU_SOURCE
// Author: Wes Kendall
// Copyright 2012 www.mpitutorial.com
// This code is provided freely with the tutorials on mpitutorial.com. Feel
// free to modify it for your own use. Any distribution of the code must
// either provide a link to www.mpitutorial.com or keep this header intact.
//
// Program that computes the average of an array of elements in parallel using
// MPI_Scatter and MPI_Gather
//
#include <stdlib.h>
#include <stdio.h>
#include "includes/dcdplugin.c"
#include "includes/priority_queue.c"
#include "includes/linked_list.c"
#include <assert.h>
// #include <omp.h>
#include "includes/mergeSort.c"
#include "includes/dynam_arr.c"
#include "includes/utility.c"

#include <time.h>
#include <mpi.h>
// #include <assert.h>


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

    // biggest and smallest value sorted by x
    int p1, p2; 
    getMinMax(Ax, &p1, &p2,fullSet, full); 
    
    // biggest and smallest value sorted by y
    int p3, p4; 
    getMinMax(Ay, &p3, &p4, fullSet, full); 

    // biggest and smallest value sorted by z
    int p5, p6; 
    getMinMax(Az, &p5, &p6,fullSet, full); 
    

    //  Find distance of each point in p array from p1 and put its square in the d1 array. 
    for (int i = 0; i < full; ++i) 
        d1[i] = getSquaredDistance(Ax[p1], Ax[fullSet[i][0]], Ay[p1], Ay[fullSet[i][0]], Az[p1], Az[fullSet[i][0]]);

    //  Find distance of each point in p array from p2 and put its square in the d2 array. 
    for (int i = 0; i < full; ++i)
        d2[i] = getSquaredDistance(Ax[p2], Ax[fullSet[i][0]], Ay[p2], Ay[fullSet[i][0]], Az[p2], Az[fullSet[i][0]]);

    //  Find distance of each point in p array from p2 and put its square in the d2 array. 
    for (int i = 0; i < full; ++i) 
        d3[i] = getSquaredDistance(Ax[p3], Ax[fullSet[i][0]], Ay[p3], Ay[fullSet[i][0]], Az[p3], Az[fullSet[i][0]]);

    //  Find distance of each point in p array from p2 and put its square in the d2 array. 
    for (int i = 0; i < full; ++i) 
        d4[i] = getSquaredDistance(Ax[p4], Ax[fullSet[i][0]], Ay[p4], Ay[fullSet[i][0]], Az[p4], Az[fullSet[i][0]]);

    //  Find distance of each point in p array from p2 and put its square in the d2 array. 
    for (int i = 0; i < full; ++i) 
        d5[i] = getSquaredDistance(Ax[p5], Ax[fullSet[i][0]], Ay[p5], Ay[fullSet[i][0]], Az[p5], Az[fullSet[i][0]]);;
    
    //  Find distance of each point in p array from p2 and put its square in the d2 array. 
    for (int i = 0; i < full; ++i) 
        d6[i] = getSquaredDistance(Ax[p6], Ax[fullSet[i][0]], Ay[p6], Ay[fullSet[i][0]], Az[p6], Az[fullSet[i][0]]);

    // using a heuristic for the sum array 
    for (int i = 0; i < full; ++i)
        sum[i] = 11*d1[i] + 101*d2[i] + 547*d3[i] + 1009*d4[i] + 5501*d5[i] + 10007*d6[i];

    // filling the index array with all the indexes
    for (int i = 0; i < full; ++i) 
        index[i] = i; 

    // merge sorting sum array and swapping corresponding indexes 
    mergeSort1(sum, index, 0, full-1); 

    // print_array_axis(sum, index, full);
    int lookAhead = 100; 
    if (k != 0) 
      lookAhead *= sqrt(k);  

    // int lookAhead = 100 * sqrt(k); 

    // TODO: priority queue with fixed size 
    float smallest = 999999999.9;
    int pA, pB; 
    for (int i = 0; i < full; ++i) {
        for (int a = 1; a <= lookAhead && a+i < full; ++a) {
            int j = a+i; 
            int pointA = fullSet[index[i]][0];
            int pointB = fullSet[index[j]][0];  
            if(fullSet[index[i]][1] != fullSet[index[j]][1]) {
                // printf("item p1: %d %d ; item p2: %d %d\n", fullSet[index[i]][0], fullSet[index[i]][1], fullSet[index[j]][0], fullSet[index[j]][1]);
                float newPoint = getSquaredDistance(Ax[pointA], Ax[pointB], Ay[pointA], Ay[pointB], Az[pointA], Az[pointB]);

                if (fullSet[index[i]][1] == 0) {
                    pushQ(pq, fullSet[index[i]][0], fullSet[index[j]][0], newPoint, k);
                } else {
                    pushQ(pq, fullSet[index[j]][0], fullSet[index[i]][0], newPoint, k);
                }
            }
        }
    }
    
    
    // printf("Smallest %d %d %lf", pA, pB, sqrtf(smallest));
    
}


int main(int argc, char** argv) {
  // if (argc != 2) {
  //   fprintf(stderr, "Usage: avg num_elements_per_proc\n");
  //   exit(1);
  // }

  int num_elements_per_proc = -1;
  // Seed the random number generator to get different results each time
  // srand(time(NULL));

  MPI_Init(NULL, NULL);

  int world_rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);
  int world_size;
  MPI_Comm_size(MPI_COMM_WORLD, &world_size);

  // DONE: calculate num_elements_per_proc using dcd->size and world_size; 
  // TODO: use scatterv instead because we need them split up nicely to all procs. 
  // nows the time to read in the input_file value 
  // DataPack ** dataPacks  = NULL;
  
    char *input_file; 
    char *output_file; 
    // getting input file name
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
    // opening the input file 
    ptr_file = fopen(input_file, "r");

    if (!ptr_file) {
        printf("NO SUCH FILE AS THE GIVEN INPUT\n"); 
        return 1;
    }

    // reading in the set values 
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
            int cS = 0; 


            dcdfile = malloc(sizeof(char*) * strlen(buf));
            //TODO: change this before you upload this dink
             while (buf[c] != '\0') {
                if (buf[c] != '\n') {
                    dcdfile[cS++] = buf[c];
                }
                c++;
             }

             dcdfile[cS-1] = '\0';

             // printf("%s", dcdfile); 
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
    // closing the input file ! 
    fclose(ptr_file);

    // using our inputs to calculate the set sizes !
    int aC = calculateSetSize(a_begins, a_ends, a_solos); 
    int bC = calculateSetSize(b_begins, b_ends, b_solos);

    int setA[aC][2]; 
    int setB[bC][2]; 

    createSet(a_begins, a_ends, a_solos,0, setA); 
    createSet(b_begins, b_ends, b_solos,1, setB); 

    // print_array2(setB, bC); 
    
    // concatenating into one full set baby! 
    int full[aC+bC][2];  

    for (int i = 0; i < aC; i++) {
        full[i][0] = setA[i][0];
        full[i][1] = setA[i][1];
    }

    for (int i = 0; i < bC; i++) {
        full[aC+i][0] = setB[i][0];
        full[aC+i][1] = setB[i][1];
    }

    int natoms = 0;
    void *raw_data = open_dcd_read(dcdfile, "dcd", &natoms);//
    if (!raw_data) {
        fprintf(stderr, "Please enter a valid name for the dcd file \n");//
        return 1;
    }
    dcdhandle *dcd = (dcdhandle *) raw_data;//

    molfile_timestep_t timestep;
    timestep.coords = (float *) malloc(3 * sizeof(float) * natoms);
    int sets = dcd->nsets; 
    // THIS IS IMPORTANT
    // THIS CONTAINS ALL PACKS 
    
    // THIS IS IMPORTANT 

    num_elements_per_proc = sets/world_size;


    // dataPacks = malloc(num_elements_per_proc * sizeof(DataPack)); 

    int startForThisProc = (world_rank) * num_elements_per_proc; 
    int endForThisProc = (world_rank+1) * num_elements_per_proc; 

    Node **listPQ = malloc(num_elements_per_proc * sizeof(Node)); 

    int leftOvers = dcd->nsets - (world_size*num_elements_per_proc); 

    int cont = 0; 
    for (int i = 0; i < dcd->nsets; i++) {
        int rc = read_next_timestep(raw_data, natoms, &timestep);
        if (rc) {
            fprintf(stderr, "error in read_next_timestep on frame %d\n", i);
            return 1;
        }
        if ( i >= startForThisProc && i < endForThisProc) {
          // printf("THIS IS HAPPENING HERE");
          calculateDistances3D(dcd->x, dcd->y, dcd->z, k-1, full, aC+bC, &(listPQ[cont]));
          cont++; 
        }
        
        if (world_rank == world_size-1) {
           if (i >= endForThisProc && i <= dcd->nsets) {
            // printf("SO THIS HAPPENED");
            calculateDistances3D(dcd->x, dcd->y, dcd->z, k-1, full, aC+bC, &(listPQ[cont]));
            cont++; 
           }
        }
    }

    MPI_Barrier(MPI_COMM_WORLD); 
    free(timestep.coords);    
    close_file_read(raw_data);

    FILE *output;
    
    for (int i = 0; i < world_size; i++) {
      MPI_Barrier(MPI_COMM_WORLD); 

      if (i == world_rank) {
        if (world_rank == 0)
          output = fopen(output_file, "w");
        else 
          output = fopen(output_file, "a+");
        // printf("Its me proc %d", world_rank);
        for (int i = 0; i < num_elements_per_proc; i++) {

          while (!isEmpty(&(listPQ[i]))) { 
            // printf("%d\n", i+(startForThisProc));

            Node *pk = peek(&(listPQ[i]));
            // printf("This is happening") ;
            fprintf(output, "\n%d, %d, %d, %f", i+startForThisProc, pk->a, pk->b, sqrtf(pk->priority)); 

            // printf("\n%d, %d, %d, %f", i+startForThisProc, pk->a, pk->b, sqrtf(pk->priority));
            pop(&(listPQ[i])); 
          } 
          // MPI_Barrier(MPI_COMM_WORLD); 
        }
        if (world_rank == world_size-1) {
          for (int i = num_elements_per_proc; i < num_elements_per_proc + leftOvers; i++) {
             while (!isEmpty(&(listPQ[i]))) { 
              // printf("%d\n", i+(startForThisProc));

              Node *pk = peek(&(listPQ[i]));
              // printf("This is happening") ;
              fprintf(output, "\n%d, %d, %d, %f", i+startForThisProc, pk->a, pk->b, sqrtf(pk->priority)); 

              // printf("\n%d, %d, %d, %f", i+startForThisProc, pk->a, pk->b, sqrtf(pk->priority));
              pop(&(listPQ[i])); 
            } 
          }
        }

        fclose(output); 
      }
      // MPI_Barrier(MPI_COMM_WORLD);  
    }
    
  MPI_Barrier(MPI_COMM_WORLD);
  MPI_Finalize();
}
