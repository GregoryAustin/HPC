/*
    Calculates the median given two lists of floats 
*/
void calculateMedian(float *A, struct LNode *R, struct LNode *G, float *median) {
    if (G == NULL && R == NULL)
        return; 

    int temp = *(int *) R->data; 
    // printf("%f\n", A[temp]); 

    // Node *tq = newNode(temp, A[temp]);

    Node *pq = newNode(temp, A[temp]); 
    int count = 1; 
    if (R != NULL) R = R->next;
    while (R != NULL){
        temp = *(int *) R->data;
        // printf("%f\n", A[temp]); 
        pushQ(&pq, temp, A[temp]); 

        // pushQ(&tq, temp, A[temp]); 

        count += 1; 
        R = R->next;
    }
    temp = *(int *)G->data;

    pushQ(&pq, temp, A[temp]); 

    // pushQ(&tq, temp, A[temp]); 

    count++; 

    if (G != NULL) G = G->next;
    while (G != NULL) {
          temp = *(int *) G->data;

          pushQ(&pq, temp, A[temp]);

          // pushQ(&tq, temp, A[temp]); 

          count += 1; 
          G = G->next;
    }

    // while (!isEmpty(&tq)) { 
    //     Node *pk = peek(&tq);
    //     printf("%d, %f\n", pk->a, pk->priority); 
    //     pop(&tq); 
    // } 

    if (count % 2 == 0) {
        // get item and (count/2) and (count/2+1) and avg them 
        
        for (int i = 0; i < (count/2) -1 ; i++) {
            pop(&pq); 
        }
        Node *pk = peek(&pq);
        float a = pk->priority; 
        pop(&pq); 
        pk = peek(&pq); 
        float b = pk->priority; 

        *median = (a+b) / 2; 
    } else {
        // count/2 - 0.5 is the middle of the list 
        float c = (float)count/2 - 0.5; 

        count = (int) c; 

        for (int i = 0; i < count; ++i) {
            pop(&pq );
        }
        Node *pk = peek(&pq);
        // printf("%d, %f\n", pk->a, pk->priority); 
        *median = pk->priority; 
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

            pushQ(pq, rInt, sqrtf(squaredDistance));

            tempG = tempG->next; 
        }

        tempR = tempR->next; 
    }
}


void bichromaticPairs3D(float *Ax, float *Ay, float *Az, int k, 
                        struct LNode* R, struct LNode* G, Node **pq) {
    if (k == 0) {

        calculateDistances3D(Ax, Ay, Az, k, R, G, pq);
        // printf ("\n\nK == 0 mate!!\nR\n");
        // printList(R, printInt); 
        // printf("\nEnd of R\n G\n");
        // printList(G, printInt); 
        // printf("\nEnd of G\n");
    } else {
        // printf("THIS IS R: \n");
        // printList(R, printInt);
        // printf("\nEND OF R\n");

        // printf("THIS IS G: \n");
        // printList(G, printInt);
        // printf("\nEND OF G\n");
        float *axis; 
        float median; 

        if (k == 1) axis = Ax; 
        else if (k == 2) axis = Ay; 
        else if (k == 3) axis = Az; 

        calculateMedian(axis, R, G, &median); // MEDIAN IS LEGIT 


        struct LNode * RL = NULL; struct LNode * RR = NULL; 
        struct LNode * GL = NULL; struct LNode * GR = NULL; 


        struct LNode * tempR = R; 
        struct LNode * tempG = G; 


        Node *Rq = NULL ;
        Node *Gq = NULL ; 
        // printf("\nFIRST LIST\n");
        // printList(R, printInt); 
        while (tempR != NULL)  {
            int tempInt = *(int *) tempR->data;

            pushQ(&Rq, tempInt, axis[tempInt]);

            tempR = tempR->next; 
        }   
        R = NULL; 
        while (!isEmpty(&Rq) ) { 
            Node *pk = peek(&Rq);
            push (&R, &(pk->a), int_size); 
                
            pop(&Rq); 
        } 
        // printf("\nSORTED LIST\n");
        // printList(R, printInt); 

        // printf("\nFIRST LIST\n");
        // printList(G, printInt); 
        while (tempG != NULL)  {
            int tempInt = *(int *) tempG->data;

            pushQ(&Gq, tempInt, axis[tempInt]);

            tempG = tempG->next; 
        }   
        G = NULL; 
        while (!isEmpty(&Gq) ) { 
            Node *pk = peek(&Gq);
            push (&G, &(pk->a), int_size); 
                
            pop(&Gq); 
        } 


        struct LNode * temp = R; 
        
        while (temp != NULL) {
            int tempInt = *(int *) temp->data;
            if (axis[tempInt] <= median) {
                // printf ("%f is smaller than %f \n", axis[tempInt], median); 
                push(&RL, &tempInt, int_size);
            } else {
                // printf ("%f is bigger than %f \n", axis[tempInt], median); 
                push(&RR, &tempInt, int_size);
            }
            temp = temp->next; 
        }
        struct LNode * RLn = NULL; struct LNode * RRn = NULL; 
        struct LNode * GLn = NULL; struct LNode * GRn = NULL; 


        struct LNode * left = RL; 
        struct LNode * right = RR;  

        while (left != NULL ) {
            int tempIntL = *(int *) left->data;
            push(&RLn, &tempIntL, int_size);
            left = left->next; 
        }

         while ( right != NULL ) {
            int tempIntR = *(int *) right->data;
            push(&RRn, &tempIntR, int_size);    
            right = right->next; 
        }

        // // set GL and GR 
        temp = G; 
        while (temp != NULL) {
            int tempInt = *(int *) temp->data;
            if (axis[tempInt] <= median) {
                push(&GL, &tempInt, int_size);
            } else {
                push(&GR, &tempInt, int_size);
            }
            temp = temp->next; 
        }

        left = GL; 
        right = GR;  
        while (left != NULL  ) {
            int tempIntL = *(int *) left->data;
            push(&GLn, &tempIntL, int_size);
            left = left->next; 
        }

        while (right != NULL ) {
            int tempIntR = *(int *) right->data;
            push(&GRn, &tempIntR, int_size);
            right = right->next; 
        }

        // printf ("\n\nCALLING RL AND GR: \n");
        if (RL != NULL && GR != NULL) bichromaticPairs3D(Ax, Ay, Az, k-1, RLn, GRn, pq);
        // printf(" DONE CALLING RL AND GR \n\n");

        // printf ("\n\nCALLING RL AND GL: ");
        if (RL != NULL && GL != NULL) bichromaticPairs3D(Ax, Ay, Az, k, RLn, GLn, pq);
        // printf(" DONE CALLING RL AND GR \n\n");
        // printf ("\n\nCALLING RR AND GR: ");
        if (RR != NULL && GR != NULL) bichromaticPairs3D(Ax, Ay, Az, k, RRn, GRn, pq);
        // printf(" DONE CALLING RL AND GR \n\n");
    }
    
}