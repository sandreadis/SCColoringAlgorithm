/* !!!IMPORTANT!!! FOLLOW THE INSTRUCTION ON LINE 102 IF YOU WANT TO RUN THE CODE ON A MTX FILE WITHOUT VALUES */

#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <string.h>
#include <omp.h>
#include <time.h>
#include <sys/time.h>
#include <stdbool.h>
#include <pthread.h>
#include "Matrix_Market/mmio.h"


#define NUM_THREADS 4

struct COO{
    int* rows;
    int* cols;
    int nodes;
    int edges;
}coo;

struct CSC_CSR{
    int* pointers;
    int* indices;
    int nodes;
    int nz;

}csc, csr;


struct COO create_coo(struct COO coo, int* I, int* J, int M, int nz);
struct CSC_CSR coo_to_csc(struct COO coo);
struct CSC_CSR coo_to_csr(struct COO coo);
void *first_trim(void* id);
void *second_trim(void* id);
void *coloring_init(void* id);
void *coloring(void* id);
void *unique_init(void* id);
void *pred_scc(void* id);

int global_scc;
bool global_check, global_check_colors, global_check_pred;
bool* used_nodes;
int* colors;
int* unique_colors;
bool* queue;
int trivial, global_color;
int iterations;

pthread_mutex_t trivial_mutex;



int main(int argc, char** argv){

    int ret_code;
    MM_typecode matcode;
    FILE *f;
    int M, N, nz;   
    int i, *I, *J;
    double *val;

    if (argc < 2)
	{
		fprintf(stderr, "Usage: %s [martix-market-filename]\n", argv[0]);
		exit(1);
	}
    else    
    { 
        if ((f = fopen(argv[1], "r")) == NULL) 
            exit(1);
    }

    if (mm_read_banner(f, &matcode) != 0)
    {
        printf("Could not process Matrix Market banner.\n");
        exit(1);
    }

    if (mm_is_complex(matcode) && mm_is_matrix(matcode) && 
            mm_is_sparse(matcode) )
    {
        printf("Sorry, this application does not support ");
        printf("Market Market type: [%s]\n", mm_typecode_to_str(matcode));
        exit(1);
    }

    /* find out size of sparse matrix .... */

    if ((ret_code = mm_read_mtx_crd_size(f, &M, &N, &nz)) !=0)
        exit(1);

    I = (int *) malloc(nz * sizeof(int));
    J = (int *) malloc(nz * sizeof(int));
    val = (double *) malloc(nz * sizeof(double));

    for (i=0; i<nz; i++)
    {
        fscanf(f, "%d %d %lg\n", &I[i], &J[i], &val[i]);
        I[i]--;  /* adjust the fscanf with fscanf(f, "%d %d\n", &I[i], &J[i]) */
        J[i]--;
    }
    free(val);

    if (f !=stdin) fclose(f);

    mm_write_banner(stdout, matcode);
    mm_write_mtx_crd_size(stdout, M, N, nz);

/******************Read matrix********************************/
    double start_time, end_time;
    start_time = omp_get_wtime();
    pthread_t threads[NUM_THREADS];
    int keep_id[NUM_THREADS];
    void *first_trim(void* id);
    
    double start, end;
    start = omp_get_wtime();
    coo = create_coo(coo, I, J, M, nz);
    csc = coo_to_csc(coo);
    csr = coo_to_csr(coo);
    free(I);
    free(J);

    trivial = 0;
    global_scc = 0;
    used_nodes = calloc(csc.nodes, sizeof(bool));
    colors = malloc(csc.nodes * sizeof(int));
    unique_colors = malloc(csc.nodes * sizeof(int));
    queue = calloc(csc.nodes, sizeof(bool));

    iterations = csc.nodes / NUM_THREADS;

    pthread_mutex_init(&trivial_mutex, NULL);

/*************************trivials***********************/
    for(int i = 0; i < NUM_THREADS; i++){
        keep_id[i] = i;
        pthread_create(&threads[i], NULL, first_trim, &keep_id[i]);
    }
    for(int i = 0; i < NUM_THREADS; i++)
        pthread_join(threads[i], NULL);
    
    for(int i = 0; i < NUM_THREADS; i++){
        keep_id[i] = i;
        pthread_create(&threads[i], NULL, second_trim, &keep_id[i]);
    }
    for(int i = 0; i < NUM_THREADS; i++)
        pthread_join(threads[i], NULL);

    free(csr.pointers);
    free(csr.indices);
    pthread_mutex_destroy(&trivial_mutex);
    global_scc = trivial;
    printf("trivial_scc: %d\n", global_scc);
/******************************************************/
    int id;
    global_check = true;
    while(global_check){
        global_check = false;
        /*****************coloring************************/
        global_check_colors = true;
        for(int i = 0; i < NUM_THREADS; i++){
            keep_id[i] = i;
            pthread_create(&threads[i], NULL, coloring_init, &keep_id[i]);
        }
        for(int i = 0; i < NUM_THREADS; i++)
            pthread_join(threads[i], NULL);

        while(global_check_colors){
            global_check_colors = false;
            for(int i = 0; i < NUM_THREADS; i++){
                keep_id[i] = i;
                pthread_create(&threads[i], NULL, coloring, &keep_id[i]);
            }
            for(int i = 0; i < NUM_THREADS; i++)
                pthread_join(threads[i], NULL);
        }
        /***************************************************/
        /*********************unique************************/
        for(int i = 0; i < NUM_THREADS; i++){
                keep_id[i] = i;
                pthread_create(&threads[i], NULL, unique_init, &keep_id[i]);
            }
        for(int i = 0; i < NUM_THREADS; i++)
            pthread_join(threads[i], NULL);
        id = 0;
        for(int i = 0; i < csc.nodes; i++){
            if(used_nodes[i] == false && colors[i] == i){
                unique_colors[id++] = i;
            }
        }
        /**************************************************/
        /*********************predecessor******************/
        for(int j = 0; j < csc.nodes; j++)
            if(unique_colors[j] != -1){
                global_color = unique_colors[j];
                queue[global_color] = true;
                used_nodes[global_color] = true;
                global_check_pred = true;
                while(global_check_pred){
                    global_check_pred = false;
                    for(int i = 0; i < NUM_THREADS; i++){
                        keep_id[i] = i;
                        pthread_create(&threads[i], NULL, pred_scc, &keep_id[i]);
                    }
                    for(int i = 0; i < NUM_THREADS; i++)
                        pthread_join(threads[i], NULL);
                }
                global_scc++;
            }
            else
                break;
                    
    /********************************************************/                
        for(int i = 0; i < csc.nodes; i++){
            if(used_nodes[i] == false){
                global_check = true;
                break;
            }
        }
    }


    free(csc.indices);
    free(csc.pointers);
    free(queue);
    free(colors);
    free(unique_colors);
    free(used_nodes);

    end_time = omp_get_wtime();
    printf("global scc: %d\ntime: %lf\n", global_scc, end_time - start_time);

    return 0;
}


void *pred_scc(void* id){

    int my_id = *(int*)id;
    int my_start, my_end, start, end, ids;
    my_start = my_id * iterations;
    my_end = (my_id + 1) * iterations;
    if(my_id = NUM_THREADS - 1) my_end = csc.nodes;
    for(int i = my_start; i < my_end; i++){
        if(queue[i] == true){
            queue[i] = false;
            start = csc.pointers[i];
            end = csc.pointers[i + 1];
            for(int j = start; j < end; j++){   
                ids = csc.indices[j];
                if(used_nodes[ids] == false && colors[ids] == global_color){
                        queue[ids] = true;
                        used_nodes[ids] = true;
                        global_check_pred = true;
                }
            }
        }
    }

}

void *unique_init(void* id){

    int my_id = *(int*)id;
    int start, end;
    start = my_id * iterations;
    end = (my_id + 1) * iterations;
    if(my_id == NUM_THREADS - 1) end = csc.nodes;
    for(int i = start; i < end; i++)
        unique_colors[i] = -1;
}

/************coloring*********/
void *coloring(void* id){

    int my_id = *(int*)id;
    int my_start, my_end, start, end, min, index;
    my_start = my_id * iterations;
    my_end = (my_id + 1) * iterations;
    if(my_id = NUM_THREADS - 1) my_end = csc.nodes;
    for(int i = my_start; i < my_end; i++){
        start = csc.pointers[i];
        end = csc.pointers[i + 1];
        min = colors[i];
        if(used_nodes[i] == false)
            for(int j = start; j < end; j++){   
                index = csc.indices[j];
                if(used_nodes[index] == false)
                    if(min > colors[index]){
                        min = colors[index];
                        colors[i] = min; 
                        global_check_colors = true;
                    }

            }
    }
}

void *coloring_init(void* id){

    int my_id = *(int*)id;
    int start, end;
    start = my_id * iterations;
    end = (my_id + 1) * iterations;
    if(my_id == NUM_THREADS - 1) end = csc.nodes;
    for(int i = start; i < end; i++)
        colors[i] = i;
}
/*****************************/

/************trivial**********/
void *first_trim(void* id){

    int my_trivials = 0;
    int my_id = *(int*) id;
    int start = my_id * iterations;
    int end = (my_id + 1) * iterations;
    if(my_id == NUM_THREADS - 1) {
        end = csc.nodes;
    }

    for(int i = start; i < end; i++){
        if(csc.pointers[i + 1] - csc.pointers[i] == 0 || csr.pointers[i + 1] - csr.pointers[i] == 0){
            used_nodes[i] = true;
            my_trivials++;
        }
    }
    pthread_mutex_lock(&trivial_mutex);
    trivial += my_trivials;
    pthread_mutex_unlock(&trivial_mutex);
}

void *second_trim(void* id){

    int my_trivials = 0;
    int my_id = *(int*) id;
    int start = my_id * iterations;
    int end = (my_id + 1) * iterations;
    if(my_id == NUM_THREADS - 1) end = csc.nodes;


    for(int j = start; j < end; j++){ 
        if(used_nodes[j] == false){
            if (csr.pointers[j + 1] - csr.pointers[j] == 1 && used_nodes[csr.indices[csr.pointers[j]]] == true){
                used_nodes[j] = true;
                my_trivials++;
            }
            else if(csc.pointers[j + 1] - csc.pointers[j] == 1 && used_nodes[csc.indices[csc.pointers[j]]] == true){
                used_nodes[j] = true;
                my_trivials++;
            }
        }
    }
    pthread_mutex_lock(&trivial_mutex);
    trivial += my_trivials;
    pthread_mutex_unlock(&trivial_mutex);
}
/*****************************/

struct CSC_CSR coo_to_csr(struct COO coo){


    struct CSC_CSR csr;
    csr.nodes = coo.nodes;
    csr.nz = coo.edges;
    csr.pointers = malloc((csr.nodes + 1) * sizeof(int));
    csr.indices = malloc(csr.nz * sizeof(int));
    
    for(int i = 0; i < csr.nodes; i++)
        csr.pointers[i] = 0;

    for(int i = 0; i < csr.nz; i++)
        csr.indices[i] = 0;



    for(int i = 0; i < csr.nz; i++)
        csr.pointers[coo.rows[i]]++;

    int temp = 0;
    int sum = 0;
    for(int i = 0; i < csr.nodes; i++){
        temp = csr.pointers[i];
        csr.pointers[i] = sum;
        sum += temp;
    }
    csr.pointers[csr.nodes] = csr.nz;

    int row = 0;
    int dest = 0;
    
    for(int i = 0; i < csr.nz; i++){
        row = coo.rows[i];
        dest = csr.pointers[row];
        csr.indices[dest] = coo.cols[i];
        csr.pointers[row]++;
    }
    
    temp = 0;
    int last = 0;
    
    for(int i = 0; i < csr.nodes + 1; i++){
        temp = csr.pointers[i];
        csr.pointers[i] = last;
        last = temp;
    }

    return csr;
}

struct CSC_CSR coo_to_csc(struct COO coo){


    struct CSC_CSR csc;
    csc.nodes = coo.nodes;
    csc.nz = coo.edges;
    csc.pointers = malloc((csc.nodes + 1) * sizeof(int));
    csc.indices = malloc(csc.nz * sizeof(int));
    
    for(int i = 0; i < csc.nodes; i++)
        csc.pointers[i] = 0;

    for(int i = 0; i < csc.nz; i++)
        csc.indices[i] = 0;



    for(int i = 0; i < csc.nz; i++)
        csc.pointers[coo.cols[i]]++;

    int temp = 0;
    int sum = 0;
    for(int i = 0; i < csc.nodes; i++){
        temp = csc.pointers[i];
        csc.pointers[i] = sum;
        sum += temp;
    }
    csc.pointers[csc.nodes] = csc.nz;

    int row = 0;
    int dest = 0;
    
    for(int i = 0; i < csc.nz; i++){
        row = coo.cols[i];
        dest = csc.pointers[row];
        csc.indices[dest] = coo.rows[i];
        csc.pointers[row]++;
    }
    
    temp = 0;
    int last = 0;
    
    for(int i = 0; i < csc.nodes + 1; i++){
        temp = csc.pointers[i];
        csc.pointers[i] = last;
        last = temp;
    }

    return csc;
} 

struct COO create_coo(struct COO coo, int* I, int* J, int M, int nz){

    coo.nodes = M;
    coo.edges = nz;
    coo.cols = J;
    coo.rows = I;
    
    return coo;
}























