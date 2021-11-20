#include <omp.h> 
#include <stdio.h> 
#include <stdlib.h> 
#include <time.h>

void crout0(double const **A, double **L, double **U, int n) {
    int i, j, k;
    double sum = 0;
    for (i = 0; i < n; i++) {
        U[i][i] = 1;
    }
    for (j = 0; j < n; j++) {
        for (i = j; i < n; i++) {
            sum = 0;
            for (k = 0; k < j; k++) {
                sum = sum + L[i][k] * U[k][j];    
            }
            L[i][j] = A[i][j] - sum;
        }
        for (i = j; i < n; i++) {
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
}

void crout1(double const **A, double **L, double **U, int n, int t) {
    #pragma omp parallel for num_threads(t)
    for (int i = 0; i < n; i++) {
        U[i][i] = 1;
    }
    for (int j = 0; j < n; j++) {
        #pragma omp parallel for num_threads(t)
        for (int i = j; i < n; i++) {
            double sum = 0;
            for (int k = 0; k < j; k++) {
                sum = sum + L[i][k] * U[k][j];    
            }
            L[i][j] = A[i][j] - sum;
        }
        #pragma omp parallel for num_threads(t)
        for (int i = j; i < n; i++) {
            double sum = 0;
            for(int k = 0; k < j; k++) {
                sum = sum + L[j][k] * U[k][i];
            }
            if (L[j][j] == 0) {                
                exit(0);
            }
            U[j][i] = (A[j][i] - sum) / L[j][j];
        }
    }
}

void crout2(double const **A, double **L, double **U, int n, int t) {
    int i, j, k;
    double sum = 0;
    omp_set_nested(1);
    omp_set_num_threads(t);
    for (i = 0; i < n; i++) {
        U[i][i] = 1;
    }
    for (j = 0; j < n; j++) {
        sum = 0;
        for (int k = 0; k < j; k++) {
            sum = sum + L[j][k] * U[k][j];    
        }
        L[j][j] = A[j][j] - sum;
        #pragma omp parallel sections
        {
        #pragma omp section
        {
        for (int i = j+1;i < n; i++) {
            double sum = 0;
            for (int k = 0; k < j; k++) {
                sum = sum + L[i][k] * U[k][j];    
            }
            L[i][j] = A[i][j] - sum;
        }
        }
        #pragma omp section
        {
        for (int i = j; i < n; i++) {
            double sum = 0;
            for(int k = 0; k < j; k++) {
                sum = sum + L[j][k] * U[k][i];
            }
            if (L[j][j] == 0) {                
                exit(0);
            }
            U[j][i] = (A[j][i] - sum) / L[j][j];
        }
        }
        }
    }
}
void crout3(double const **A, double **L, double **U, int n, int t) {
    omp_set_nested(1);
    omp_set_num_threads(t);
    int i, j;
    #pragma omp parallel for num_threads(t)
    for (i = 0; i < n; i++) {
        U[i][i] = 1;
    }
    for (j = 0; j < n; j++) {
        double sum = 0;
        //#pragma omp parallel for num_threads(t)
        for (int k = 0; k < j; k++) {
            sum = sum + L[j][k] * U[k][j];    
        }
        L[j][j] = A[j][j] - sum;
        #pragma omp parallel sections
        {
        #pragma omp section
        {
        #pragma omp parallel for
        for (int i = j+1; i < n; i++) {
            double sum = 0;
            for (int k = 0; k < j; k++) {
                sum = sum + L[i][k] * U[k][j];    
            }
            L[i][j] = A[i][j] - sum;
        }
        //L[j...n][j]      //U[0...j][j]

        //U[j][j..n-1]...L[j][j]
        }
        #pragma omp section
        {
        #pragma omp parallel for
        for (int i = j; i < n; i++) {
            double sum = 0;
            for(int k = 0; k < j; k++) {
                sum = sum + L[j][k] * U[k][i];
            }
            if (L[j][j] == 0) {                
                exit(0);
            }
            
        U[j][i] = (A[j][i] - sum) / L[j][j];
        }
        }
        }
    }




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
    //omp_set_nested(1);
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
    
    if(strategy==0){

        crout0(A,L,U,n);
        write_output(output_L,L,n);
        write_output(output_U,U,n);
    }
    else if(strategy==1){
        crout1(A,L,U,n,t);
        write_output(output_L,L,n);
        write_output(output_U,U,n);
    }
    else if(strategy==2){
        crout2(A,L,U,n,t);
        write_output(output_L,L,n);
        write_output(output_U,U,n);
    }
    else if(strategy==3){
        crout3(A,L,U,n,t);
        write_output(output_L,L,n);
        write_output(output_U,U,n);                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                       
    }
    else{

    }
    //printf("%lf",check(L,U,A,n));
    free(L);
    free(U);
    free(A);
    fclose(input);
    return 0;
}


