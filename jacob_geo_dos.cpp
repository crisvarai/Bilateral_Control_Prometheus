#include <iostream>
#include <cmath>
#include <iomanip>
#include <chrono>

using namespace std;
using namespace std::chrono;

//----------------- Declaración de variables del Jacobiano Geométrico -----------------------//

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

double sig_32, sig_31;
double sig_30, sig_29, sig_28, sig_27, sig_26, sig_25, sig_24, sig_23, sig_22, sig_21;
double sig_20, sig_19, sig_18, sig_17, sig_16, sig_15, sig_14, sig_13, sig_12, sig_11;
double sig_10, sig_9,  sig_8,  sig_7,  sig_6,  sig_5,  sig_4,  sig_3,  sig_2,  sig_1;

double Jv_11, Jv_21, Jv_31, Jv_12, Jv_22, Jv_32, Jv_13, Jv_23, Jv_33, Jv_14, Jv_24, Jv_34, Jv_15, Jv_25, Jv_35, Jv_16, Jv_26, Jv_36;
double Jw_41, Jw_51, Jw_61, Jw_42, Jw_52, Jw_62, Jw_43, Jw_53, Jw_63, Jw_44, Jw_54, Jw_64, Jw_45, Jw_55, Jw_65, Jw_46, Jw_56, Jw_66;

void jacob_geo(double q1, double q2, double q3, double q4, double q5, double q6, double a2, double d3, double a3, double a4, double d5, double d7, double k) {
    //------------- Variables sigma para reducir la expresión del Jacobiano -------------//
    sig_32 = sin(k)*sin(q2) + cos(k)*cos(q1)*cos(q2);
    sig_31 = cos(k)*sin(q2) - cos(q1)*cos(q2)*sin(k);
    sig_30 = cos(q1)*sin(q3) + cos(q2)*cos(q3)*sin(q1);
    sig_29 = cos(k)*cos(q1)*sin(q2);
    sig_28 = cos(q2)*sin(k);
    sig_27 = cos(q3)*sig_32 - cos(k)*sin(q1)*sin(q3);
    sig_26 = cos(k)*cos(q3)*sin(q1);
    sig_25 = cos(k)*cos(q2) + cos(q1)*sin(k)*sin(q2);
    sig_24 = cos(q3)*sig_31 + sin(k)*sin(q1)*sin(q3);
    sig_23 = sin(q4)*sig_30 + cos(q4)*sin(q1)*sin(q2);
    sig_22 = cos(q4)*sig_30 - sin(q1)*sin(q2)*sin(q4);
    sig_21 = cos(q1)*cos(q3) - cos(q2)*sin(q1)*sin(q3);
    sig_20 = sin(q4)*sig_27 - cos(q4)*(sig_28 - sig_29);
    sig_19 = cos(q4)*sig_27 + sin(q4)*(sig_28 - sig_29);
    sig_18 = sin(q3)*sig_32 + sig_26;
    sig_17 = sin(q4)*sig_24 - cos(q4)*sig_25;
    sig_16 = d7*cos(q6)*sig_23;
    sig_15 = d7*sin(q6)*(cos(q5)*sig_22 + sin(q5)*sig_21);
    sig_14 = a4*cos(q4)*sig_30;
    sig_13 = a4*sin(q1)*sin(q2)*sin(q4);
    sig_12 = cos(q4)*sig_24 + sin(q4)*sig_25;
    sig_11 = sin(q3)*sig_31 - cos(q3)*sin(k)*sin(q1);
    sig_10 = d7*cos(q6)*sig_20;
    sig_9  = sin(q5)*sig_22 - cos(q5)*sig_21;
    sig_8  = cos(q5)*sig_18 + sin(q5)*sig_19;
    sig_7  = d7*cos(q6)*sig_17;
    sig_6  = d5*sig_23 + sig_16 + sig_15;
    sig_5  = d5*sig_23 + sig_16 + sig_15 + sig_14 - sig_13;
    sig_4  = cos(q5)*sig_11 + sin(q5)*sig_12;
    sig_3  = d7*sin(q6)*(sin(q5)*sig_18 - cos(q5)*sig_19);
    sig_2  = d7*sin(q6)*(sin(q5)*sig_11 - cos(q5)*sig_12);
    sig_1  = d5*sig_23 + a2*cos(q2)*sin(q1) + a3*cos(q1)*sin(q3) + d3*sin(q1)*sin(q2) + sig_16 + sig_15 + sig_14 + a3*cos(q2)*cos(q3)*sin(q1) - sig_13;

    //------------------- Ecuaciones para indexar en la triz del Jacobiano --------------------//

    // Sección que representa la velocindad lineal
    Jv_11 = -cos(k)*sig_1;
    Jv_21 = a2*cos(q1)*cos(q2) + d3*cos(q1)*sin(q2) - a3*sin(q1)*sin(q3) + a3*cos(q1)*cos(q2)*cos(q3) + d5*cos(q1)*cos(q4)*sin(q2) - a4*cos(q1)*sin(q2)*sin(q4) - a4*cos(q4)*sin(q1)*sin(q3) - d5*sin(q1)*sin(q3)*sin(q4) + a4*cos(q1)*cos(q2)*cos(q3)*cos(q4) + d5*cos(q1)*cos(q2)*cos(q3)*sin(q4) + d7*cos(q1)*cos(q4)*cos(q6)*sin(q2) - d7*cos(q6)*sin(q1)*sin(q3)*sin(q4) - d7*cos(q3)*sin(q1)*sin(q5)*sin(q6) - d7*cos(q1)*cos(q2)*sin(q3)*sin(q5)*sin(q6) - d7*cos(q1)*cos(q5)*sin(q2)*sin(q4)*sin(q6) - d7*cos(q4)*cos(q5)*sin(q1)*sin(q3)*sin(q6) + d7*cos(q1)*cos(q2)*cos(q3)*cos(q6)*sin(q4) + d7*cos(q1)*cos(q2)*cos(q3)*cos(q4)*cos(q5)*sin(q6);
    Jv_31 = -sin(k)*sig_1;

    Jv_12 = a2*cos(q2)*sin(k) + d3*sin(k)*sin(q2) + d3*cos(k)*cos(q1)*cos(q2) - a2*cos(k)*cos(q1)*sin(q2) + a3*cos(q2)*cos(q3)*sin(k) + d5*cos(q4)*sin(k)*sin(q2) - a4*sin(k)*sin(q2)*sin(q4) + d5*cos(k)*cos(q1)*cos(q2)*cos(q4) - a3*cos(k)*cos(q1)*cos(q3)*sin(q2) - a4*cos(k)*cos(q1)*cos(q2)*sin(q4) + a4*cos(q2)*cos(q3)*cos(q4)*sin(k) + d5*cos(q2)*cos(q3)*sin(k)*sin(q4) + d7*cos(q4)*cos(q6)*sin(k)*sin(q2) - d5*cos(k)*cos(q1)*cos(q3)*sin(q2)*sin(q4) + d7*cos(q2)*cos(q3)*cos(q6)*sin(k)*sin(q4) - d7*cos(q2)*sin(k)*sin(q3)*sin(q5)*sin(q6) - d7*cos(q5)*sin(k)*sin(q2)*sin(q4)*sin(q6) + d7*cos(k)*cos(q1)*cos(q2)*cos(q4)*cos(q6) - a4*cos(k)*cos(q1)*cos(q3)*cos(q4)*sin(q2) - d7*cos(k)*cos(q1)*cos(q3)*cos(q6)*sin(q2)*sin(q4) - d7*cos(k)*cos(q1)*cos(q2)*cos(q5)*sin(q4)*sin(q6) + d7*cos(q2)*cos(q3)*cos(q4)*cos(q5)*sin(k)*sin(q6) + d7*cos(k)*cos(q1)*sin(q2)*sin(q3)*sin(q5)*sin(q6) - d7*cos(k)*cos(q1)*cos(q3)*cos(q4)*cos(q5)*sin(q2)*sin(q6);
    Jv_22 = -sin(q1)*(a2*sin(q2) - d3*cos(q2) - d5*cos(q2)*cos(q4) + a3*cos(q3)*sin(q2) + a4*cos(q2)*sin(q4) - d7*cos(q2)*cos(q4)*cos(q6) + a4*cos(q3)*cos(q4)*sin(q2) + d5*cos(q3)*sin(q2)*sin(q4) + d7*cos(q3)*cos(q6)*sin(q2)*sin(q4) + d7*cos(q2)*cos(q5)*sin(q4)*sin(q6) - d7*sin(q2)*sin(q3)*sin(q5)*sin(q6) + d7*cos(q3)*cos(q4)*cos(q5)*sin(q2)*sin(q6));
    Jv_32 = d3*cos(q1)*cos(q2)*sin(k) - d3*cos(k)*sin(q2) - a3*cos(k)*cos(q2)*cos(q3) - a2*cos(k)*cos(q2) - d5*cos(k)*cos(q4)*sin(q2) - a2*cos(q1)*sin(k)*sin(q2) + a4*cos(k)*sin(q2)*sin(q4) - a4*cos(k)*cos(q2)*cos(q3)*cos(q4) + d5*cos(q1)*cos(q2)*cos(q4)*sin(k) - d5*cos(k)*cos(q2)*cos(q3)*sin(q4) - d7*cos(k)*cos(q4)*cos(q6)*sin(q2) - a3*cos(q1)*cos(q3)*sin(k)*sin(q2) - a4*cos(q1)*cos(q2)*sin(k)*sin(q4) - a4*cos(q1)*cos(q3)*cos(q4)*sin(k)*sin(q2) - d5*cos(q1)*cos(q3)*sin(k)*sin(q2)*sin(q4) + d7*cos(k)*cos(q2)*sin(q3)*sin(q5)*sin(q6) + d7*cos(k)*cos(q5)*sin(q2)*sin(q4)*sin(q6) + d7*cos(q1)*cos(q2)*cos(q4)*cos(q6)*sin(k) - d7*cos(k)*cos(q2)*cos(q3)*cos(q6)*sin(q4) + d7*cos(q1)*sin(k)*sin(q2)*sin(q3)*sin(q5)*sin(q6) - d7*cos(k)*cos(q2)*cos(q3)*cos(q4)*cos(q5)*sin(q6) - d7*cos(q1)*cos(q3)*cos(q6)*sin(k)*sin(q2)*sin(q4) - d7*cos(q1)*cos(q2)*cos(q5)*sin(k)*sin(q4)*sin(q6) - d7*cos(q1)*cos(q3)*cos(q4)*cos(q5)*sin(k)*sin(q2)*sin(q6);

    Jv_13 = d7*cos(k)*sin(q1)*sin(q3)*sin(q5)*sin(q6) - a3*sin(k)*sin(q2)*sin(q3) - a3*cos(k)*cos(q1)*cos(q2)*sin(q3) - a4*cos(k)*cos(q3)*cos(q4)*sin(q1) - d5*cos(k)*cos(q3)*sin(q1)*sin(q4) - a4*cos(q4)*sin(k)*sin(q2)*sin(q3) - d5*sin(k)*sin(q2)*sin(q3)*sin(q4) - d5*cos(k)*cos(q1)*cos(q2)*sin(q3)*sin(q4) - d7*cos(k)*cos(q3)*cos(q6)*sin(q1)*sin(q4) - a3*cos(k)*cos(q3)*sin(q1) - d7*cos(q6)*sin(k)*sin(q2)*sin(q3)*sin(q4) - d7*cos(q3)*sin(k)*sin(q2)*sin(q5)*sin(q6) - a4*cos(k)*cos(q1)*cos(q2)*cos(q4)*sin(q3) - d7*cos(k)*cos(q1)*cos(q2)*cos(q6)*sin(q3)*sin(q4) - d7*cos(k)*cos(q1)*cos(q2)*cos(q3)*sin(q5)*sin(q6) - d7*cos(k)*cos(q3)*cos(q4)*cos(q5)*sin(q1)*sin(q6) - d7*cos(q4)*cos(q5)*sin(k)*sin(q2)*sin(q3)*sin(q6) - d7*cos(k)*cos(q1)*cos(q2)*cos(q4)*cos(q5)*sin(q3)*sin(q6);
    Jv_23 = a3*cos(q1)*cos(q3) + a4*cos(q1)*cos(q3)*cos(q4) + d5*cos(q1)*cos(q3)*sin(q4) - a3*cos(q2)*sin(q1)*sin(q3) + d7*cos(q1)*cos(q3)*cos(q6)*sin(q4) - a4*cos(q2)*cos(q4)*sin(q1)*sin(q3) - d5*cos(q2)*sin(q1)*sin(q3)*sin(q4) - d7*cos(q1)*sin(q3)*sin(q5)*sin(q6) - d7*cos(q2)*cos(q6)*sin(q1)*sin(q3)*sin(q4) - d7*cos(q2)*cos(q3)*sin(q1)*sin(q5)*sin(q6) + d7*cos(q1)*cos(q3)*cos(q4)*cos(q5)*sin(q6) - d7*cos(q2)*cos(q4)*cos(q5)*sin(q1)*sin(q3)*sin(q6);
    Jv_33 = a3*cos(k)*sin(q2)*sin(q3) - a3*cos(q3)*sin(k)*sin(q1) - a3*cos(q1)*cos(q2)*sin(k)*sin(q3) - a4*cos(q3)*cos(q4)*sin(k)*sin(q1) + a4*cos(k)*cos(q4)*sin(q2)*sin(q3) - d5*cos(q3)*sin(k)*sin(q1)*sin(q4) + d5*cos(k)*sin(q2)*sin(q3)*sin(q4) - a4*cos(q1)*cos(q2)*cos(q4)*sin(k)*sin(q3) - d5*cos(q1)*cos(q2)*sin(k)*sin(q3)*sin(q4) - d7*cos(q3)*cos(q6)*sin(k)*sin(q1)*sin(q4) + d7*cos(k)*cos(q6)*sin(q2)*sin(q3)*sin(q4) + d7*cos(k)*cos(q3)*sin(q2)*sin(q5)*sin(q6) + d7*sin(k)*sin(q1)*sin(q3)*sin(q5)*sin(q6) - d7*cos(q1)*cos(q2)*cos(q6)*sin(k)*sin(q3)*sin(q4) - d7*cos(q1)*cos(q2)*cos(q3)*sin(k)*sin(q5)*sin(q6) - d7*cos(q3)*cos(q4)*cos(q5)*sin(k)*sin(q1)*sin(q6) + d7*cos(k)*cos(q4)*cos(q5)*sin(q2)*sin(q3)*sin(q6) - d7*cos(q1)*cos(q2)*cos(q4)*cos(q5)*sin(k)*sin(q3)*sin(q6);

    Jv_14 = -sig_11*sig_5 - sig_21*(d5*sig_17 + a4*cos(q4)*sig_24 - sig_2 + sig_7 + a4*sin(q4)*sig_25);
    Jv_24 = d5*cos(q1)*cos(q4)*sin(q3) - a4*cos(q4)*sin(q1)*sin(q2) - a4*cos(q1)*sin(q3)*sin(q4) - d5*sin(q1)*sin(q2)*sin(q4) + d5*cos(q2)*cos(q3)*cos(q4)*sin(q1) + d7*cos(q1)*cos(q4)*cos(q6)*sin(q3) - a4*cos(q2)*cos(q3)*sin(q1)*sin(q4) - d7*cos(q6)*sin(q1)*sin(q2)*sin(q4) - d7*cos(q4)*cos(q5)*sin(q1)*sin(q2)*sin(q6) - d7*cos(q1)*cos(q5)*sin(q3)*sin(q4)*sin(q6) + d7*cos(q2)*cos(q3)*cos(q4)*cos(q6)*sin(q1) - d7*cos(q2)*cos(q3)*cos(q5)*sin(q1)*sin(q4)*sin(q6);
    Jv_34 = -sig_18*sig_5 - sig_21*(d5*sig_20 - sig_3 + a4*cos(q4)*sig_27 + sig_10 + a4*sin(q4)*(sig_28 - sig_29));

    Jv_15 = sig_17*sig_6 - sig_23*(d5*sig_17 - sig_2 + sig_7);
    Jv_25 = -d7*sin(q6)*(cos(q2)*cos(q5)*sin(q1)*sin(q3) - cos(q1)*cos(q3)*cos(q5) + cos(q1)*cos(q4)*sin(q3)*sin(q5) - sin(q1)*sin(q2)*sin(q4)*sin(q5) + cos(q2)*cos(q3)*cos(q4)*sin(q1)*sin(q5));
    Jv_35 = sig_20*sig_6 - sig_23*(d5*sig_20 - sig_3 + sig_10);

    Jv_16 = -(sig_16 + sig_15)*sig_4 - sig_9*(sig_2 - sig_7);
    Jv_26 = (sig_2 - sig_7)*sig_8 - (sig_3 - sig_10)*sig_4;
    Jv_36 = -sig_8*(sig_16 + sig_15) - sig_9*(sig_3 - sig_10);

    // Sección que representa la velocindad angular
    Jw_41 = -sin(k);
    Jw_51 = 0;
    Jw_61 = cos(k);
    
    Jw_42 = - cos(k)*sin(q1);
    Jw_52 = cos(q1);
    Jw_62 = -sin(k)*sin(q1);

    Jw_43 = sig_29 - sig_28;
    Jw_53 = sin(q1)*sin(q2);
    Jw_63 = sig_25;

    Jw_44 = -sin(q3)*sig_32 - sig_26;
    Jw_54 = sig_21;
    Jw_64 = sig_11;

    Jw_45 = sig_20;
    Jw_55 = sig_23;
    Jw_65 = cos(q4)*sig_25 - sin(q4)*sig_24;

    Jw_46 = -cos(q5)*sig_18 - sin(q5)*sig_19;
    Jw_56 = cos(q5)*sig_21 - sin(q5)*sig_22;
    Jw_66 = sig_4;

    double J[6][6] = {{Jv_11,Jv_12,Jv_13,Jv_14,Jv_15,Jv_16},
                      {Jv_21,Jv_22,Jv_23,Jv_24,Jv_25,Jv_26},
                      {Jv_31,Jv_32,Jv_33,Jv_34,Jv_35,Jv_36},
                      {Jw_41,Jw_42,Jw_43,Jw_44,Jw_45,Jw_46},
                      {Jw_51,Jw_52,Jw_53,Jw_54,Jw_55,Jw_56},
                      {Jw_61,Jw_62,Jw_63,Jw_64,Jw_65,Jw_66}};
        
    // cout << "Jacobiano :" << endl;
    // for(int x=0;x<6;x++)  
    // {
    //   for(int y=0;y<6;y++) 
    //     {
    //       cout << setw(10) << setprecision(4) << fixed << J[x][y]; // setw() = dimensión de la salida, fixed() y setprecision() = decimales
    //     }
    // cout<<endl;
    // }
}

int main()
{
    auto passes = 0;
    auto tStart = high_resolution_clock::now();

    while(true)
    {
        jacob_geo(q1,q2,q3,q4,q5,q6,a2,d3,a3,a4,d5,d7,k);
        auto tEnd = high_resolution_clock::now();
        duration<float> duration = tEnd - tStart;
        passes++;
        if (duration.count() >= 1)
        {
            cout << setprecision(2) << fixed << "Duration :" << duration.count() << endl;
            cout << setprecision(6) << fixed << "Passes :" << passes << endl;
            break;
        }
    }
    return 0;
}
