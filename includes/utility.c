// #include "dynam_arr.c"
#include <stdlib.h>
#include <stdio.h>


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

/*
    Creates a array[] set 
*/
void createSet(Array begins, Array ends, Array solos, int setID, int list[][2]) {    

    int count = 0; 
    int size = 0; 

    for (int i = 0; i < begins.used; ++i){ 
        for (int j = 0; j <= (ends.array[i] - begins.array[i]); j++) {
            list[size+j][0] = begins.array[i]+j; 
            list[size+j][1] = setID;
        }
        size += ends.array[i] - begins.array[i] + 1; 
    }

    for (int i = 0; i < solos.used; ++i) {
        list[size+i][0] = solos.array[i];
        list[size+i][1] = setID; 
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

int calculateSetSize(Array begins, Array ends, Array solos) {
    int count = 0; 

    for (int i = 0; i < begins.used; ++i) 
        count +=  ends.array[i] - begins.array[i] + 1;

    for (int i = 0; i < solos.used; ++i) 
        count += 1; 

    return count; 
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