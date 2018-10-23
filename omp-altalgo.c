#define _GNU_SOURCE
#include <stdlib.h>
#include <stdio.h>
#include "includes/dcdplugin.c"
#include "priority_queue.c"
#include "linked_list.c"
#include <assert.h>
#include <time.h>
#include <omp.h>
#include "mergeSort.c"


typedef int bool;
#define true 1
#define false 0

#define NUM_THREADS 8

unsigned int_size = sizeof(int); 

void print_array(int *array, int length)
{
    for (int i = 0; i < length; i++) { 
        printf("%d ", array[i]);

     }
     printf("\n");
}

void print_array2(int array[][2], int length)
{
    for (int i = 0; i < length; i++) { 
        printf("%d %d, ", array[i][0], array[i][1]);

     }
     printf("\n");
}

void print_array_axis(float * axis, int *array, int length)
{
    for (int i = 0; i < length; i++) { 
        printf("%f \n", axis[array[i]]);

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

bool checkDifferentSets(int a, int b, int list[], int size) {
    bool aBool = false; 
    bool bBool = false; 

    for (int i = 0; i < size; i++) {
        if (a == list[i]) 
            aBool = true; 
        if (b == list[i])
            bBool = true; 

        if (bBool == true && aBool == true) {
            // printf("\n%d same set %d\n", a, b);
            return false; 
        }
    }

    if (aBool == true && bBool == false || aBool == false && bBool == true) 
        return true; 
    else 
        return false; 
    
}

/*
    Creates a LinkedList set given starting index and ending idnex  
*/
void createSet(int set_begin, int set_end, int setID, int list[][2]) {    

    int i; 
    for (i = 0; i <= (set_end - set_begin); i++) {
        list[i][0] = set_begin+i; 
        list[i][1] = setID; 

    }
    
}

void *array_concat(const void *a, size_t an, const void *b, size_t bn, size_t s)
{
    char *p = malloc(s * (an + bn));
    memcpy(p, a, an*s);
    memcpy(p + an*s, b, bn*s);
    return p;
}

float getSquaredDistance(float dxA, float dxB, float dyA, float dyB, float dzA, float dzB) {
    float dx = dxA - dxB;
    float dy = dyA - dyB;
    float dz = dzA - dzB;

    float dx2 = dx * dx;
    float dy2 = dy * dy;
    float dz2 = dz * dz;

    float squaredDistance = dx2 + dy2 + dz2;

    return squaredDistance;
}



void calculateDistances3D(float *Ax, float *Ay, float *Az, int k, int fullSet[][2], int full, Node **pq) {
    float d1[full], d2[full], d3[full], d4[full], d5[full], d6[full], sum[full];
    int index[full];

    // biggest and smallest value sorted by x
    mergeSort2(Ax, fullSet, 0, full-1);
    int p1 = fullSet[0][0]; 
    int p2 = fullSet[full-1][0];

    // biggest and smallest value sorted by y
    mergeSort2(Ay, fullSet, 0, full-1);
    int p3 = fullSet[0][0];
    int p4 = fullSet[full-1][0];

    // biggest and smallest value sorted by z
    mergeSort2(Az, fullSet, 0, full-1);
    int p5 = fullSet[0][0];
    int p6 = fullSet[full-1][0];

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
    

    // TODO: priority queue with fixed size 
    float smallest = 999999999.9;
    int pA, pB; 
    for (int i = 0; i < full; ++i) {
        for (int j = i+1; j <= i+100 && j < full; ++j) {
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

  	int setA[aC][2]; 
    int setB[bC][2]; 

    
    // printf("About to create sets%d", fSize);;
    createSet(a_begin, a_end, 0, setA); 
    createSet(b_begin, b_end, 1, setB); 

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
        

        // print_array_axis(dcd->x, full, aC + bC);
        // float *Ax, float *Ay, float *Az, int k, int fullSet[], int listA[], int aC, int bC
        calculateDistances3D(dcd->x, dcd->y, dcd->z, k, full, aC+bC, &pq);  

        int count = 0; 
        while (!isEmpty(&pq) && count < 3) { 
            Node *pk = peek(&pq);
            printf("%d, %d, %d, %f\n", i, pk->a, pk->b, sqrtf(pk->priority)); 
            pop(&pq); 
            count++;
        } 

        // cleaning up! 
        while (!isEmpty(&pq)) 
            pop(&pq); 
    
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