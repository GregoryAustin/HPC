#define _GNU_SOURCE
#include <stdlib.h>
#include <stdio.h>
#include "includes/dcdplugin.c"
#include "priorityqueue.c"
#include "linked_list.c"
#include <assert.h>
#include <time.h>
#include <math.h>

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

unsigned float_size = sizeof(float); 
unsigned int_size = sizeof(int); 

typedef int bool;
#define true 1
#define false 0

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
void createSet(int set_begin, int set_end, struct LNode** head_ref) {
    int i; 
    
    for (i = set_end; i >= set_begin; --i) {
        // printf("PUSHING"); 
        push(head_ref, &i, int_size); 
    }
}

void sort(float *A, struct LNode *R, struct LNode *G, struct LNode **RG, int * count) {
    Node *pq = NULL; 
    int temp; 
    *count = 0; 
    // putting all of R into pq 
    while (R != NULL) {
        temp = *(int *) R->data;
        pushQ(&pq, temp, A[temp]); 
        R = R->next;
    }

    // putting all of G into pq 
    while (G != NULL) {
          temp = *(int *) G->data;
          pushQ(&pq, temp, A[temp]);
          G = G->next;
    }

    // putting the sorted list into RG 
    while (!isEmpty(&pq)) { 
        Node *pk = peek(&pq);
        push(RG, &(pk->a), int_size);
        pop(&pq); 
        *count += 1; 
    } 
}


void flipLinkedList(struct LNode **lst) {
    struct LNode * temp = *lst; 
    struct LNode * tmpNode = NULL;

    while (temp != NULL ) {
        int tInt = *(int *) temp->data;
        push(&tmpNode, &tInt, int_size);
        temp = temp->next; 
    }

    *lst = tmpNode; 
}


void print_array(int *array, int length)
{
    for (int i = 0; i < length; i++) { 
        printf("%d ", array[i]);

     }
     printf("\n");
}


bool checkDifferentSets(int a, int b, struct LNode* R) {
    bool aBool = false; 
    bool bBool = false; 
    struct LNode * temp = R; 

    while (temp != NULL) {
        int tempInt = *(int *) temp->data;
        if (a == tempInt) 
            aBool = true; 
        if (b == tempInt)
            bBool = true; 

        if (bBool == true && aBool == true) {
            // printf("\n%d same set %d\n", a, b);
            return false; 
        }

        temp = temp->next; 
    }

    if (aBool == true && bBool == false || aBool == false && bBool == true) {
        // printf("\n%d diff set %d\n", a, b);
        return true; 
    }
    else {
        // printf("\n%d same set %d\n", a, b);
        return false; 
    }
}

float calculate (float *Ax, float *Ay, float *Az, int a, int b) {
    float dx = Ax[a] - Ax[b];
    float dy = Ay[a] - Ay[b];
    float dz = Az[a] - Az[b];

    float dx2 = dx * dx;
    float dy2 = dy * dy;
    float dz2 = dz * dz;

    float squaredDistance = dx2 + dy2 + dz2;
    return sqrtf(squaredDistance);
}

float closest(float *Ax, float *Ay, float *Az, int listx[], int size,
    struct LNode* R) { 
    
    float closestDist = 999999.0;
     // printf ("\n hi\n");
    printf("\n size is %d ", size); 
    printf ("\nChecking this biz out..  \n");
    print_array(listx, size); 

    for (int i = 0; i < size; i++) {
        for (int j = i + 1; j < size; j++) {
            // TODO: check you're comparing set with other set
            float dist; 
            printf("%d %d \n",i, listx[i]);
            printf("%d %d \n",i, listx[j]);
            if (checkDifferentSets(listx[i], listx[j], R)) 
                dist = calculate(Ax, Ay, Az, listx[i], listx[j]);
            else break; 

            if (dist < closestDist) 
                closestDist = dist;
        }
    }
    printf ("\nI made it through \n");
    return closestDist; 
}

float midlaneClosest(float *Ax, float *Ay, float *Az, 
    int midlane[], int length, float smallest,
    struct LNode* R) {

    printf ("HOWZIT!"); 
    // TODO: check you're comparing set with other set 

    float newMin = smallest; 

    for (int i = 0; i < length; i++) {
        for (int j = i + 1; (j < length) && ((Ay[midlane[j]] - Ay[midlane[i]]) < newMin); j++) {
            float min; 
            if (checkDifferentSets(midlane[i], midlane[j], R))
                min = calculate(Ax, Ay, Az, midlane[i], midlane[j]);
            else break; 

            if (min < newMin) 
                newMin = min; 
        }
    }

    return newMin; 
}



float divideAndConquer(float *Ax, float *Ay, float *Az, 
    int listByX [],  int listByY [], int count, 
    struct LNode* R) {
    // printf("Check me baby!\n");
    if (count <= 3) {
        // brute force as per usual 
        printf("Im here baby!\n");


        
        return closest(Ax, Ay, Az, listByX, count, R);
    }



    printf("\nAnd the game begins... \n");
    printf("COUNT IS %d\n", count); 
    int intermediate = ceil(count / 2.0);//(int)ceil(count / 2.0); 
    printf("INTERMEDIATE IS %d\n", intermediate);
    print_array(listByX, count); 
    int midpoint = listByX[intermediate-1];
    printf("MIDPOINT IS %d\n", midpoint);

    int PL[intermediate+1];
    int PR[count - intermediate - 1]; 

    int QL[intermediate+1];
    int QR[count - intermediate - 1]; 

    for (int a = 0; a < count; a++) {

        if (a <= intermediate) {
            PL[a] = listByX[a]; 
            QL[a] = listByY[a];
        }
        else {
            PR[a] = listByX[a]; 
            QR[a] = listByY[a];
        }
     
    }   
    printf("Check me baby!\n");

    float length1 = divideAndConquer(Ax, Ay, Az, PL, QL, intermediate, R); 
    float length2 = divideAndConquer(Ax, Ay, Az, PR, QR, count-intermediate, R); 
    printf("After that recursion  baby!\n");
    float smallest;

    if (length1 < length2)
        smallest = length1;
    else 
        smallest = length2; 

    int midlane[count]; 
    int midlength = 0; 
    printf("Still going baby!\n");
    printf("%d is the Count\n", count);
    printf("%d is the midlength\n", midlength); 
    for (int a = 0; a < count; a++) {
        printf ("%d ITEM start \n", a); 
        printf ("%d list by y \n", listByY[a]);
        printf ("%d midpoint \n", midpoint);
        printf ("%f list by y \n", Ax[listByY[a]]);
        printf ("%f midpoint \n", Ax[midpoint]);
        printf ("%f smallest \n", smallest);

        if (abs(Ax[listByY[a]] - Ax[midpoint]) < smallest) {
            printf("here"); 
            midlane[midlength] = listByY[a];
            midlength++; 
        }
        printf ("%d ITEM end \n", a);
    }
    printf("Still at it baby!\n");
    float midl = midlaneClosest(Ax, Ay, Az, midlane, midlength, smallest, R);
    if (smallest < midl) 
        return smallest;
    else 
        return midl;

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

    // TODO: adapt this for ranges, and itemized lists 
    struct LNode *setA = NULL; 
    struct LNode *setB = NULL; 

    createSet(a_begin, a_end, &setA); // THESE ARE 100s 
    createSet(b_begin, b_end, &setB); // THESE ARE 100s 

    /***************************************************************/
    /**************** Done Reading Input File **********************/
    /***************************************************************/
    
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
        Node *pq = NULL; //newNode(1,999999); 
        // bichromaticPairs3D(dcd->x, dcd->y, dcd->z,  3, setB, setA, &pq);
        
        // TODO:
        // sort setA union setB according to X coordinate 
        struct LNode *listByX = NULL; 
        struct LNode *listByY = NULL;

        // int count = 0; 
        sort(dcd->x, setA, setB, &listByX, &count);
        

        // // sort setA union setB according to Y coordinate 
        sort(dcd->y, setA, setB, &listByY, &count);

        int listX[count];
        int listY[count];   
        printList(listByX, printInt); 
        int c = 0;
        while (listByX != NULL) {
            listX[c] = *(int *) listByX->data;
            listY[c] = *(int *) listByY->data;

            listByX = listByX->next; 
            listByY = listByY->next; 
            c++; 
        }   

        
        printf("\n\n\n");
        print_array(listX, count);

        printf("%d closest points: %f\n", i, divideAndConquer(dcd->x, dcd->y, dcd->z, listX, listY, count, setA)); 

        // int c = 0; 
        // while (!isEmpty(&pq) && c < 3) { 
        //     Node *pk = peek(&pq);
        //     printf("%d, %d, %f\n", i, pk->a, pk->priority); 
        //     pop(&pq); 
        //     c++;
        // } 
        
        // Node *pq = newNode(1,1,999999); 
        // // printf("HERE");
        // calculateDistances3D(dcd->x, dcd->y, dcd->z, 3, a_begin, a_end, b_begin, b_end, &pq);  


        // int count = 0; 
        // while (!isEmpty(&pq) && count < 3) { 
        //     Node *pk = peek(&pq);
        //     printf("%d, %d, %d, %f\n", i, pk->a, pk->b, pk->priority); 
        //     pop(&pq); 
        //     count++;
        // } 

        // // cleaning up! 
        // while (!isEmpty(&pq)) 
        //     pop(&pq); 
    
        // printf("\n");


        // printf("Timestep %d\n", i);
        // printf("i: x    y      z\n");
        
        if (i >= 0) break;
    }
    free(timestep.coords);    
    close_file_read(raw_data);
    
    return 0;

}