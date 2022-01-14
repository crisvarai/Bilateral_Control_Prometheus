#include <iostream>
#include <cmath>    // To use "pi" and "trigonometric" functions 
#include <iomanip>  // To "setw" and "setprecision" functions

using namespace std;

//----------------- Declaración de variables para calcular la Traza -----------------------//
double matrix[4][4]={{2,3,5,1},{7,4,5,1},{8,5,4,1},{1,1,1,1}};
double traza = 0;
const int N = 4;

//----------------- Declaración de funciones para calcular la Traza -----------------------//
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
            cout << setw(8) << setprecision(4) << fixed << arr[i][j];
        }
        cout << endl;
    }
}

double trace(double** mat) {
    for(int i=0;i<N;i++){
        for(int j=0;j<N;j++) {
            if(i==j)                        // condition for trace of matrix
                traza=traza+mat[i][j];
        }
    }
    return traza;
}

//------------------------------- INICIO DEL PROGRAMA PRINCIPAL -----------------------------//
int main()
{
    double **arr;
    arr = getArray(matrix);
    cout << "matrix: " << endl;
    printArray(arr);
    
    double result = trace(arr);
    cout << "trace: " << result << endl;
    return 0;
}