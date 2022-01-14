// Return two dimensional array from function in C/C++
// Method: Using Dynamic Array

#include <iostream>

using namespace std;

const int N = 3;

void printArray(int** arr) {
    for(int i=0; i<N; i++) {
        for(int j=0; j<N; j++) {
            cout << arr[i][j];
        }
        cout << endl;
    }
}

int** getArray() {
    int** arr = new int*[N];
    for(int i=0; i<N; i++) {
        arr[i] = new int[N];
        for(int j=0; j<N; j++) {
            arr[i][j] = i+j;
        }
    }
    return arr;
}

int main() 
{
    int **arr;
    arr = getArray();
    printArray(arr);
    return 0;
}
  