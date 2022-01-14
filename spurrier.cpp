#include <iostream>
#include <cmath>     // To use "pi" and "trigonometric" functions 
#include <iomanip>   // To "setw" and "setprecision" functions
#include <complex>   // To manipulate complex numbers "imag"
#include <algorithm> // To find max element of a vector

using namespace std;

//----------------- Declaración de variables para Quaternion (Spurrier) -----------------------//
const int H = 4;
const int N = 3;
double x = 0, y = 0, z = 0;
double w = 0;

double T[4][4] = {{0.4226, 0.0000, -0.9063, -595.0771}, 
                  {0.0000, 1.0000,  0.0000,    0.0000},
                  {0.9063, 0.0000,  0.4226,  382.9389},
                  {0,      0,       0,         1     }};

// double T[4][4] = {{-1.0000, -0.0000, -0.0000, -224.7959}, 
//                   {-0.0000,  1.0000,  0.0000,  -37.2900}, 
//                   { 0.0000,  0.0000, -1.0000, -439.7823}, 
//                   { 0,       0,       0,        1     }};

// double T[4][4] = {{-0.3712, 0.0000, -0.9285, -679.70}, 
//                   { 0.0000, 1.0000,  0.0000,  -22.24},
//                   { 0.9285, 0.0000, -0.3712,  -32.52},
//                   { 0,      0,       0,         1  }};

//----------------- Declaración de funciones para Quaternion (Spurrier) -----------------------//
double** getArray(double mat[H][H]) {
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
    double traza = 0;
    for(int i=0;i<N;i++){
        for(int j=0;j<N;j++) {
            if(i==j)                        // condition for trace of matrix
                traza=traza+mat[i][j];
        }
    }
    return traza;
}

double* getVector(double vec[H]) {
    double* arr = new double[H];
    for(int i=0;i<H;i++) {
        arr[i] = vec[i];
    }
    return arr;
}

void printVector(double* arr) {
    for(int i=0; i<H; i++) {
        cout << setw(8) << setprecision(4) << fixed << arr[i];
    } 
    cout << endl;
}

double* spurrier(double mat[H][H]) {
    double **R = getArray(T);
    // cout << "matrix: " <<endl;
    // printArray(R);
    
    double w = trace(R);
    // cout << "trace: " << w << endl;

    double img = imag(w);
    if(img>0) { 
        w = 0; 
    };
    // cout << "img: " << img << endl;
    // cout << "trace: " << w << endl;

    double vec[4] = {w, R[0][0], R[1][1], R[2][2]};
    // cout << "vec: " << endl;
    // for(int i=0;i<4;i++) {
    //     cout << setw(8) << setprecision(4) << fixed << vec[i];
    // } cout << endl;

    double *maxim = max_element(begin(vec),end(vec));
    // cout << "max: " << *maxim << endl;

    int pos = distance(begin(vec),max_element(begin(vec),end(vec)));
    // cout << "pos: " << pos << endl;

    if(pos == 0) {
        w = sqrt(trace(R)+1)/2;
        x = (R[2][1]-R[1][2])/(4*w);
        y = (R[0][2]-R[2][0])/(4*w);
        z = (R[1][0]-R[0][1])/(4*w);
    }

    if(pos == 1) {
        x = sqrt((R[0][0]/2)+((1-trace(R))/4));
        w = (R[2][1]-R[1][2])/(4*x);
        y = (R[1][0]+R[0][1])/(4*x);
        z = (R[2][0]+R[0][2])/(4*x);
    }

    if(pos == 2) {
        y = sqrt((R[1][1]/2)+((1-trace(R))/4));
        w = (R[0][2]-R[2][0])/(4*y);
        z = (R[2][1]+R[1][2])/(4*y);
        x = (R[0][1]+R[1][0])/(4*y);
    }

    if(pos == 3 ) {
        z = sqrt((R[2][2]/2)+((1-trace(R))/4));
        w = (R[1][0]-R[0][1])/(4*z);
        x = (R[0][2]+R[2][0])/(4*z);
        y = (R[1][2]+R[2][1])/(4*z);
    }

    double quaternion[4] = {w,x,y,z};
    getVector(quaternion);
}

//------------------------------- INICIO DEL PROGRAMA PRINCIPAL -----------------------------//
int main()
{
    double *result;
    result = spurrier(T);
    printVector(result);
    return 0;
}