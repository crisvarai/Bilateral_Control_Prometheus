#include <iostream>
#include <string>
#include <vector>
#include <sstream>

#include <iomanip>   // To "setw" and "setprecision" functions

using namespace std;

string data = "123.45, 213.76, 234.16";

double* parse(string s) {
    size_t pos = 0;
    vector< double > vd;
    double d = 0.0;

    // convert ',' to ' '
    while (pos < s.size())
        if ((pos = s.find_first_of (',',pos)) != string::npos)
            s[pos] = ' ';

    stringstream ss(s);
    while (ss >> d)
        vd.push_back(d);

    double* arr = new double[vd.size()];
    for(int i=0; i<vd.size(); i++) {
        arr[i] = vd[i];
    }
    return arr;
}



int main () {
    double *res = parse(data);
    cout << "parse: " << endl;
    for(int i=0; i<3; i++) {
        cout << setw(8) << setprecision(4) << fixed << res[i] << "\n";
    } 
    cout << endl;
    return 0;
}