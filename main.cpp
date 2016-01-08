#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <fstream>
#include <math.h>

using namespace std;

/*------------------Declaração das funções-------------------------*/
//void mostraDados(int n_temp, int k_temp, double* y_temp, double** X_temp);
//void mostraMatriz(int fila, int coluna,  double** X_temp);
//void mostraVetor(int fila, double* b);

void produtoAb(int filaA, int colunaA, double** A, double* b, double*&C); // Produto matriz vetor : c = Ab
void LU(int n, double* b_temp, double** A_temp, double *&beta); // Resolve o sistema de equações LU beta = y
void residualAdjustado(int n, int k, double *e, double** X, double *&r);
void geraAb(int n, int k, double** X, double** W, double* y, double**&A, double*&b); // Gera a matriz A e o vetor b

#define eps 1.0e-10  // Tolerância

int main(){
    /*------------------Leitura do arquivo de dados-------------------------*/
    ifstream file;
    file.open("instances/prob_2.txt"); // mudanças para o git

    if (file.fail( )){
        printf(" Input file %s opening failed.\n ", "/instances");
        system(" pause ");
        exit(1);
    }
    
    int n, k, i, j, t;

    /*----------------- Leitura da instancia------------------*/
    file >> n;
    file >> k;
    printf("n =%d \n", n);
    printf("k =%d \n", k);
    double y[n];
    double *py = &y[0]; 
    double X[n][k];
    double **pX;
    if (file.is_open()){
        for (i = 0; i < k; i++)
            for (j = 0; j < n; j++){
                file >> X[j][i];
            }

        for (t = 0; t < n; t++)
            file >> y[t];
    }
    file.close();
    
    /* -----------------Declaração de varaveis-------------*/
    double W[n][n];  // Matriz de ponderação
    double **pW;
    double A[k+1][k+1]; 
    double **pA;   
    double b[k+1];
    double *pb = &b[0];
    
    double *p[n];
    double *ww[n];
    for(i = 0; i < n; i++){
        p[i] = X[i];
        ww[i] = W[i];
    }
    pX = p;
    pW = ww;
    
    double *aa[k];
    for(i = 0; i < k+1; i++)
        aa[i] = A[i];
    pA = aa;
    
    double beta0[k];
    double *pbeta0 = &beta0[0];
    double Xbeta[n];
    double *pXbeta = &Xbeta[0];     
    double beta1[k];
    double *pbeta1 = &beta1[0];
    
    /* -----------------Calcula o valor de beta0----------------------*/ 
    for(i = 0; i < n; i++)
        for(j = 0; j < n; j++)
            if(i == j)
                W[i][j] = 1;
            else
                W[i][j] = 0;
    
    geraAb(n, k, pX, pW, y, pA, pb);
    LU(k+1, pb, pA, pbeta0);  // Calcula beta0   
    //mostraVetor(k,beta0);
    
    /* -----------------Processo iterativo----------------------------*/ 
    double e[n], r[n], w[n];
    double *pr = &r[0];
    
    int Itermax = 100;
        
    for(i = 1; i < Itermax; i++){
        
        produtoAb(n, k, pX, pbeta0, pXbeta);
        
        for(j = 0; j < n; j++) 
            e[j] = y[j] - Xbeta[j];
        
        residualAdjustado(n, k, e, pX, pr);
        
        /* Função de ponderação bisquare*/
        for(t = 0; t < n; t++){
            r[t] = r[t] / 4.685;
            if( fabs(r[t]) < 1 )
                w[t] = pow((1-pow(r[t],2)),2);
            else
                w[t] = 0;
      
            W[t][t] = w[t];
        } 
        
        geraAb(n, k, pX, pW, y, pA, pb);
        LU(k+1, pb, pA, pbeta1);  // Calcula beta1
        //mostraVetor(k,beta1);
        /*-------------- Criterio de parada-------------------- */
        int band = 1;
        
        for(t = 0; t < k; t++){
            if(fabs(beta1[t]-beta0[t]) > eps){
                band = 0;
            }
        }
        
        if(band == 1){
            printf(" i = %d \n", i);
            double sum = 0;
            printf("A solução é beta : \n");
            for(t = 0; t < k; t++)
                printf(" beta[%d] = %.10f \n", t, beta1[t]);
            
            break;
        }
        else{
            for(t = 0; t < k; t++)
                beta0[t] = beta1[t];
        }
    } 
    
    
    return 0;
}
