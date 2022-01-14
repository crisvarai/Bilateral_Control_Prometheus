#include <iostream>
#include <cmath>    // To use "pi" and "trigonometric" functions 
#include <iomanip>  // To "setw" and "setprecision" functions

using namespace std;

//----------------- Declaración de variables para Cinemática Directa -----------------------//
double Cq, Sq, Calpha, Salpha;
double **T0B, **TBH, **T1, **T1_NEW, **T2, **T3, **T4, **T5, **T6, **TWM, **T6_NEW, **TEE;
const int N = 4;

// Declaración de parámetros DH 
double q1 = ( 0 )*M_PI/180;
double q2 = ( 0 )*M_PI/180;
double q3 = ( 0 )*M_PI/180;
double q4 = ( 0 )*M_PI/180;
double q5 = ( 0 )*M_PI/180;
double q6 = ( 0 )*M_PI/180;

double d1 = 183.73;     // dist. del pecho al hombro (q2)
double a2 = 58.28;
double d3 = 244.4;
double a3 = 18.76;
double a4 = 18.53;
double d5 = 123.03;
double d7 = 150.0;      // 150 mm is dist. from q6 to ee; 71 mm, from q7.

double k = 65*M_PI/180; // k equivale a alpha en DH

//----------------- Declaración de funciones para Cinemática Directa -----------------------//
double** getArray(double mat[N][N]) {
    // Etapa que permite regresar un arreglo de dos dimensiones (2d array)
    double** arr = new double*[N];
    for(int i=0; i<N; i++) {
        arr[i] = new double[N];
        for(int j=0; j<N; j++) {
            arr[i][j] = mat[i][j];
        }
    }
    return arr;
}

void printArray(double** arr) {
    for(int i=0; i<N; i++) {
        for(int j=0; j<N; j++) {
            cout << setw(10) << setprecision(4) << fixed << arr[i][j];
        }
        cout << endl;
    }
}

double** DH(double q, double d, double a, double alpha) {
    Cq = cos(q);
    Sq = sin(q);
    Calpha = cos(alpha);
    Salpha = sin(alpha);
    double Tzrx[N][N] = {{Cq,-Calpha*Sq,Salpha*Sq,a*Cq}, {Sq,Calpha*Cq,-Cq*Salpha,a*Sq}, {0,Salpha,Calpha,d}, {0,0,0,1}};
    getArray(Tzrx);
}

double** multiply(double** mat1, double** mat2)
{
    double res[N][N];
    for (int i = 0; i < N; i++) {
        for (int j = 0; j < N; j++) {
            res[i][j] = 0;
            for (int k = 0; k < N; k++) {
                res[i][j] += mat1[i][k] * mat2[k][j];
            }
        }
    }
    getArray(res);
}

double** FK() {
    T0B = DH(M_PI/2,0,0,-k);
    TBH = DH(-M_PI/2,0,0,0);
    T1 = DH(q1,d1,0,-M_PI/2);
    T1_NEW = multiply(multiply(T0B,TBH),T1);
    T2 = DH(q2,0,a2,M_PI/2);
    T3 = DH(q3,d3,a3,-M_PI/2);
    T4 = DH(q4,0,a4,M_PI/2);
    T5 = DH(q5,d5,0,-M_PI/2);
    T6 = DH(q6,0,0,M_PI/2);
    TWM = DH(0,d7,0,0);
    T6_NEW = multiply(T6,TWM);
    TEE = multiply(multiply(multiply(multiply(multiply(T1_NEW,T2),T3),T4),T5),T6_NEW); 
    return TEE;
}

//------------------------------- INICIO DEL PROGRAMA PRINCIPAL -----------------------------//
int main()
{
    double **result;
    result = FK();
    cout << "FK:" << endl;
    printArray(result);
    return 0;
}
