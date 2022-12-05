/* !!!IMPORTANT!!! Follow the instruction on line 95 to run an mtx file without values column */

#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <string.h>
#include <omp.h>
#include <time.h>
#include <sys/time.h>
#include "Matrix_Market/mmio.h"
#include <stdbool.h>



/*********************My structs**************************/
struct COO{
    int* rows;
    int* cols;
    int nodes;
    int edges;
};

struct CSC_CSR{
    int* pointers;
    int* indices;
    int nodes;
    int nz;

};
/**********************************************************/

/*********************My functions***********************/

int trivial_scc(struct CSC_CSR csc, struct CSC_CSR csr, bool* used_nodes);
int* coloring(struct CSC_CSR csc, bool* used_nodes);
int* unique(int* colors, bool* used_nodes, struct CSC_CSR csc);
void pred_scc(int* colors, struct CSC_CSR csc, bool* used_nodes, int index);
struct CSC_CSR coo_to_csc(struct COO coo);
struct COO create_coo(struct COO coo, int* I, int* J, int M, int nz);
struct CSC_CSR coo_to_csr(struct COO coo);


/********************************************************/

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


    /* reseve memory for matrices */

    I = (int *) malloc(nz * sizeof(int));
    J = (int *) malloc(nz * sizeof(int));
    val = (double *) malloc(nz * sizeof(double));

    for (i=0; i<nz; i++)
    {
        fscanf(f, "%d %d %lg\n", &I[i], &J[i], &val[i]);
        I[i]--;  /* change the scanf to scanf(f, "%d %d\n", &I[i], &J[i]) */
        J[i]--;
    }
    free(val);

    if (f !=stdin) fclose(f);


    mm_write_banner(stdout, matcode);
    mm_write_mtx_crd_size(stdout, M, N, nz);

/******************Read matrix********************************/

    double start, end;
    start = omp_get_wtime();
    struct COO coo = create_coo(coo, I, J, M, nz);
    struct CSC_CSR csc = coo_to_csc(coo);
    struct CSC_CSR csr = coo_to_csr(coo);
    free(I);
    free(J);

/***********************************************************/

    int global_scc, used_nodes_num = 0;
    bool* used_nodes = calloc(csc.nodes, sizeof(bool));
    int trivial = trivial_scc(csc, csr, used_nodes);
    global_scc = trivial;
    
    printf("trivial_scc: %d\n", global_scc);

    int *colors, *unique_colors;

    bool done = false;
    while(!done){
        done = true;
        colors = coloring(csc, used_nodes);
        unique_colors = unique(colors, used_nodes, csc);

        for(int i = 0; i < csc.nodes; i++){
            if(unique_colors[i] != -1){   
                pred_scc(colors, csc, used_nodes, unique_colors[i]);
                global_scc++;  
            }
            else
                break;
        }

        for(int i = 0; i < csc.nodes; i++){
            if(used_nodes[i] == false){
                done = false;
                break;
            }
        }
        free(colors);
        free(unique_colors);       
    }


    printf("global_scc: %d\n", global_scc);




    free(used_nodes);
    free(csc.pointers);
    free(csc.indices);
    free(csr.pointers);         
    free(csr.indices);
    end = omp_get_wtime();
    printf("time: %lf\n", end - start);

    return 0;
}












void pred_scc(int* colors, struct CSC_CSR csc, bool* used_nodes, int index){

    bool check = true;
    int id, start, end;
    bool* queue = calloc(csc.nodes, sizeof(bool));
    queue[index] = 1;
    used_nodes[index] = true;


    while(check){
        check = false;
        for(int i = 0; i < csc.nodes; i++){ 
            if(queue[i] == 1){
                queue[i] = 0;
                start = csc.pointers[i];
                end = csc.pointers[i + 1];
                for(int j = start; j < end; j++){
                    id = csc.indices[j];
                    if(colors[index] == colors[id] && used_nodes[id] == false){
                        queue[id] = 1;
                        check = true;
                        used_nodes[id] = true;
                    }
                }
            }
        }
    }
}

int* unique(int* colors, bool* used_nodes, struct CSC_CSR csc){

    int* unique_colors = malloc(csc.nodes * sizeof(int));
    int id = 0;

    for(int i = 0; i < csc.nodes; i++)
        unique_colors[i] = -1;

    for(int i = 0; i < csc.nodes; i++){ 
        if(used_nodes[i] == false && colors[i] == i){
            unique_colors[id] = i;
            id++;
        }
    }

    return unique_colors;
}

int* coloring(struct CSC_CSR csc, bool* used_nodes){

    bool check = true;
    int start, end, index, min;
    int* colors = malloc(csc.nodes * sizeof(int));
    for(int i = 0; i < csc.nodes; i++){ 
        if(used_nodes[i] == false)
            colors[i] = i;
        else colors[i] = -1;
    }

    while (check){
        check = false;
        for(int i = 0; i < csc.nodes; i++){
            start = csc.pointers[i];
            end = csc.pointers[i + 1];
            min = colors[i];
            if(used_nodes[i] == false){
                for(int j = start; j < end; j++){
                    index = csc.indices[j];
                    if(used_nodes[index] == false){
                        if(min > colors[index]){
                            min = colors[index];
                            colors[i] = min;
                            check = true;
                        }
                    }
                }
            }
        }
    }

    return colors;
}

int trivial_scc(struct CSC_CSR csc, struct CSC_CSR csr, bool* used_nodes){

    int trivial = 0;

    for(int j = 0; j < csc.nodes; j++){
        if(csc.pointers[j + 1] - csc.pointers[j] == 0 || csr.pointers[j + 1] - csr.pointers[j] == 0){
                used_nodes[j] = true;
                trivial = trivial + 1; 
        }
    }
    for(int j = 0; j < csc.nodes; j++){ 
        if(used_nodes[j] == false){
            if (csr.pointers[j + 1] - csr.pointers[j] == 1 && used_nodes[csr.indices[csr.pointers[j]]] == true){
                used_nodes[j] = true;
                trivial = trivial + 1; 
            }
            else if(csc.pointers[j + 1] - csc.pointers[j] == 1 && used_nodes[csc.indices[csc.pointers[j]]] == true){
                used_nodes[j] = true;
                trivial = trivial + 1;
            }
        }
    }

    return trivial;
}

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










