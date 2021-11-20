#include <stdio.h> 
#include <stdlib.h> 
#include <time.h>
#include <math.h>
#include "mpi.h"

void crout4(double const **A, double **L, double **U, int n, int p) {
    int my_rank,comm_sz;
    int i, j, k;
    /*double **local_U = (double **) malloc(n * sizeof(double*));
    double **local_L = (double **) malloc(n * sizeof(double*));*/
    MPI_Init(NULL,NULL);
    MPI_Comm_rank(MPI_COMM_WORLD,&my_rank);
    MPI_Comm_size(MPI_COMM_WORLD,&comm_sz);
    for(int i=0;i<n;i++){
        if(i%comm_sz==my_rank){
                U[i][i]=1;
        } 
    }
    for(int i=0;i<n;i++){
        MPI_Bcast(&U[i][i],1,MPI_DOUBLE,i%comm_sz,MPI_COMM_WORLD);
    }
    double sum = 0;
    for (j = 0; j < n; j++) {
        for (i = j; i < n; i++) {
            if(i%comm_sz==my_rank){
                sum = 0;
                for (k = 0; k < j; k++) {
                    sum = sum + L[i][k] * U[k][j];    
                }
                L[i][j] = A[i][j] - sum;
            }    
            
            
        }
        for(i=j;i<n;i++){
            MPI_Bcast(&L[i][j],1,MPI_DOUBLE,i%comm_sz,MPI_COMM_WORLD);
        }
        for (i = j; i < n; i++) {
            if(i%comm_sz==my_rank){
                sum = 0;
                for(k = 0; k < j; k++) {
                    sum = sum + L[j][k] * U[k][i];
                }
                if (L[j][j] == 0) {                
                    exit(0);
                }
                U[j][i] = (A[j][i] - sum) / L[j][j];
            }
            
        }
        for (i = j; i < n; i++) { 
            MPI_Bcast(&U[j][i],1,MPI_DOUBLE,i%comm_sz,MPI_COMM_WORLD);
        }
        
    }
    MPI_Finalize();



}

void write_output(char fname[], double **arr, int n ){
    FILE *f = fopen(fname, "w");
    for( int i = 0; i < n; i++){
        for(int j = 0; j < n; j++){
            fprintf(f, "%0.12f ", arr[i][j]);
        }
        fprintf(f, "\n");
    }
    fclose(f);
}

int main(int argc,char* argv[]){
    int n = atoi(argv[1]);
    char* filename = argv[2];
    int t = atoi(argv[3]);
    int strategy = atoi(argv[4]);
    FILE* input = fopen(filename,"r");
    double **A= (double **) malloc(n * sizeof(double*));
    double **L= (double **) malloc(n * sizeof(double*));
    double **U= (double **) malloc(n * sizeof(double*));
    for (int i=0; i<n; i++) 
    {
        A[i] = (double*) malloc(n* sizeof(double));
        L[i] = (double*) malloc(n* sizeof(double));
        U[i] = (double*) malloc(n* sizeof(double));

    }
    for(int i=0;i<n;i++){
        for(int j=0;j<n;j++){
            fscanf(input,"%lf",&A[i][j]);
            L[i][j]=0;
            U[i][j]=0;
        }
    }


    char output_L[100] = "output_L_";
    strcat(output_L,argv[4]);
    strcat(output_L,"_");
    strcat(output_L,argv[3]);
    strcat(output_L,".txt");
    char output_U[100] = "output_U_";
    strcat(output_U,argv[4]);
    strcat(output_U,"_");
    strcat(output_U,argv[3]);
    strcat(output_U,".txt");                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                            
    if(strategy==4){
        crout4(A,L,U,n,t);
        write_output(output_L,L,n);
        write_output(output_U,U,n);
    }
    else{

    }
    free(L);
    free(U);
    free(A);
    fclose(input);
    return 0;
}

