// Return two dimensional array from function in C/C++
// Method: Using Dynamic Array

#include <iostream>
#include <iomanip>

using namespace std;

const int N = 4;
double matrix[4][4] = {{1,2,3,4}, {5,6,7,8}, {9,10,11,12}, {13,14,15,16}};

void printArray(double** arr) {
    for(int i=0; i<N; i++) {
        for(int j=0; j<N; j++) {
            cout << setw(9) << setprecision(4) << fixed << arr[i][j];
        }
        cout << endl;
    }
}

double** getArray(double mat[N][N]) {
    double** arr = new double*[N];
    for(int i=0; i<N; i++) {
        arr[i] = new double[N];
        for(int j=0; j<N; j++) {
            arr[i][j] = mat[i][j];
        }
    }
    return arr;
}

int main() 
{
    double **arr;
    arr = getArray(matrix);
    printArray(arr);
    return 0;
}