#include <iostream>
#include <cmath>
#include <iomanip>
using namespace std;

//------------------------- Declaración de parámetros DH ----------------------------//
double q1 = ( -60 )*M_PI/180;
double q2 = ( 30 )*M_PI/180;
double q3 = ( -45 )*M_PI/180;
double q4 = ( 60 )*M_PI/180;
double q5 = ( -30 )*M_PI/180;
double q6 = ( 45 )*M_PI/180;

double d1 = 183.73;  // dist. del pecho al hombro (q2)
double a2 = 58.28;
double d3 = 244.4;
double a3 = 18.76;
double a4 = 18.53;
double d5 = 123.03;
double d7 = 150.0;     // 150 mm is dist. from q6 to ee; 71 mm, from q7.

double k = 65*M_PI/180;

double sig_66, sig_65, sig_64, sig_63, sig_62, sig_61;
double sig_60, sig_59, sig_58, sig_57, sig_56, sig_55, sig_54, sig_53, sig_52, sig_51;
double sig_50, sig_49, sig_48, sig_47, sig_46, sig_45, sig_44, sig_43, sig_42, sig_41;
double sig_40, sig_39, sig_38, sig_37, sig_36, sig_35, sig_34, sig_33, sig_32, sig_31;
double sig_30, sig_29, sig_28, sig_27, sig_26, sig_25, sig_24, sig_23, sig_22, sig_21;
double sig_20, sig_19, sig_18, sig_17, sig_16, sig_15, sig_14, sig_13, sig_12, sig_11;
double sig_10, sig_9,  sig_8,  sig_7,  sig_6,  sig_5,  sig_4,  sig_3,  sig_2,  sig_1;

double Jv_11, Jv_21, Jv_31, Jv_12, Jv_22, Jv_32, Jv_13, Jv_23, Jv_33, Jv_14, Jv_24, Jv_34, Jv_15, Jv_25, Jv_35, Jv_16, Jv_26, Jv_36;
double Jw_41, Jw_51, Jw_61, Jw_42, Jw_52, Jw_62, Jw_43, Jw_53, Jw_63, Jw_44, Jw_54, Jw_64, Jw_45, Jw_55, Jw_65, Jw_46, Jw_56, Jw_66;

void jacob_geo(double q1, double q2, double q3, double q4, double q5, double q6, double a2, double d3, double a3, double a4, double d5, double d7, double k) {

    //------------- Variables sigma para reducir la expresión del Jacobiano -------------//
    sig_66 = sin(k)*sin(q2) + cos(k)*cos(q1)*cos(q2);
    sig_65 = cos(k)*sin(q2) - cos(q1)*cos(q2)*sin(k);
    sig_64 = sin(q1)*sin(q3) - cos(q1)*cos(q2)*cos(q3);
    sig_63 = cos(k)*cos(q1)*sin(q3) + cos(k)*cos(q2)*cos(q3)*sin(q1);
    sig_62 = cos(q1)*sin(k)*sin(q3) + cos(q2)*cos(q3)*sin(k)*sin(q1);
    sig_61 = cos(q2)*sin(k) - cos(k)*cos(q1)*sin(q2);
    sig_60 = cos(k)*cos(q2) + cos(q1)*sin(k)*sin(q2);
    sig_59 = cos(q1)*sin(q3) + cos(q2)*cos(q3)*sin(q1);
    sig_58 = cos(q3)*sig_66 - cos(k)*sin(q1)*sin(q3);
    sig_57 = cos(q3)*sig_65 + sin(k)*sin(q1)*sin(q3);
    sig_56 = sin(q3)*sig_65 - cos(q3)*sin(k)*sin(q1);
    sig_55 = sin(q3)*sig_66 + cos(k)*cos(q3)*sin(q1);
    sig_54 = cos(q1)*cos(q3) - cos(q2)*sin(q1)*sin(q3);
    sig_53 = cos(q3)*sin(q1) + cos(q1)*cos(q2)*sin(q3);
    sig_52 = cos(q4)*sig_64 + cos(q1)*sin(q2)*sin(q4);
    sig_51 = cos(k)*cos(q1)*cos(q3) - cos(k)*cos(q2)*sin(q1)*sin(q3);
    sig_50 = cos(q4)*sig_63 - cos(k)*sin(q1)*sin(q2)*sin(q4);
    sig_49 = cos(q1)*cos(q3)*sin(k) - cos(q2)*sin(k)*sin(q1)*sin(q3);
    sig_48 = cos(q4)*sig_62 - sin(k)*sin(q1)*sin(q2)*sin(q4);
    sig_47 = sin(q4)*sig_66 - cos(q3)*cos(q4)*sig_61;
    sig_46 = sin(q4)*sig_65 - cos(q3)*cos(q4)*sig_60;
    sig_45 = cos(q4)*sig_59 - sin(q1)*sin(q2)*sin(q4);
    sig_44 = cos(q4)*sig_58 + sin(q4)*sig_61;
    sig_43 = cos(q4)*sig_57 + sin(q4)*sig_60;
    sig_42 = sin(q4)*sig_57 - cos(q4)*sig_60;
    sig_41 = sin(q4)*sig_58 - cos(q4)*sig_61;
    sig_40 = sin(q4)*sig_59 + cos(q4)*sin(q1)*sin(q2);
    sig_39 = sin(q5)*sig_57 + cos(q4)*cos(q5)*sig_56;
    sig_38 = sin(q5)*sig_58 + cos(q4)*cos(q5)*sig_55;
    sig_37 = cos(q2)*sin(q1)*sin(q4) + cos(q3)*cos(q4)*sin(q1)*sin(q2);
    sig_36 = sin(q5)*sig_59 - cos(q4)*cos(q5)*sig_54;
    sig_35 = cos(q5)*sig_52 + sin(q5)*sig_53;
    sig_34 = sin(q4)*sig_64 - cos(q1)*cos(q4)*sin(q2);
    sig_33 = sin(q5)*sig_51 + cos(q5)*sig_50;
    sig_32 = sin(q4)*sig_63 + cos(k)*cos(q4)*sin(q1)*sin(q2);
    sig_31 = sin(q5)*sig_49 + cos(q5)*sig_48;
    sig_30 = sin(q4)*sig_62 + cos(q4)*sin(k)*sin(q1)*sin(q2);
    sig_29 = cos(q5)*sig_47 + sin(q3)*sin(q5)*sig_61;
    sig_28 = cos(q4)*sig_66 + cos(q3)*sin(q4)*sig_61;
    sig_27 = cos(q5)*sig_46 + sin(q3)*sin(q5)*sig_60;
    sig_26 = cos(q4)*sig_65 + cos(q3)*sin(q4)*sig_60;
    sig_25 = cos(q5)*sig_45 + sin(q5)*sig_54;
    sig_24 = sin(q5)*sig_55 - cos(q5)*sig_44;
    sig_23 = sin(q5)*sig_56 - cos(q5)*sig_43;
    sig_22 = cos(q6)*sig_43 - cos(q5)*sin(q6)*sig_42;
    sig_21 = cos(q6)*sig_44 - cos(q5)*sin(q6)*sig_41;
    sig_20 = cos(q6)*sig_45 - cos(q5)*sin(q6)*sig_40;
    sig_19 = sin(q6)*sig_39 + cos(q6)*sin(q4)*sig_56;
    sig_18 = sin(q6)*sig_38 + cos(q6)*sin(q4)*sig_55;
    sig_17 = cos(q5)*sig_37 - sin(q1)*sin(q2)*sin(q3)*sin(q5);
    sig_16 = sin(q6)*sig_36 - cos(q6)*sin(q4)*sig_54;
    sig_15 = cos(q6)*sig_34 + sin(q6)*sig_35;
    sig_14 = sin(q6)*sig_33 + cos(q6)*sig_32;
    sig_13 = cos(q6)*sig_30 + sin(q6)*sig_31;
    sig_12 = cos(q6)*sig_28 - sin(q6)*sig_29;
    sig_11 = cos(q6)*sig_26 - sin(q6)*sig_27;
    sig_10 = sin(q5)*sig_45 - cos(q5)*sig_54;
    sig_9  = cos(q6)*sig_40 + sin(q6)*sig_25;
    sig_8  = cos(q5)*sig_55 + sin(q5)*sig_44;
    sig_7  = cos(q5)*sig_56 + sin(q5)*sig_43;
    sig_6  = sin(q6)*sig_40 - cos(q6)*sig_25;
    sig_5  = sin(q6)*sig_24 - cos(q6)*sig_41;
    sig_4  = sin(q6)*sig_23 - cos(q6)*sig_42;
    sig_3  = sin(q6)*sig_41 + cos(q6)*sig_24;
    sig_2  = cos(q6)*sig_23 + sin(q6)*sig_42;
    sig_1  = cos(q2)*cos(q4)*sin(q1) - cos(q3)*sin(q1)*sin(q2)*sin(q4);
    
    //------------------- Ecuaciones para indexar en la matriz del Jacobiano --------------------//

    // Sección que representa la velocindad lineal
    Jv_11 = a4*cos(k)*sin(q1)*sin(q2)*sin(q4) - d7*sig_14 - a4*cos(q4)*sig_63 - a2*cos(k)*cos(q2)*sin(q1) - a3*cos(k)*cos(q1)*sin(q3) - d3*cos(k)*sin(q1)*sin(q2) - a3*cos(k)*cos(q2)*cos(q3)*sin(q1) - d5*sig_32;
    Jv_21 = a2*cos(q1)*cos(q2) - d7*sig_15 - d5*sig_34 + d3*cos(q1)*sin(q2) - a3*sin(q1)*sin(q3) - a4*cos(q4)*sig_64 + a3*cos(q1)*cos(q2)*cos(q3) - a4*cos(q1)*sin(q2)*sin(q4);
    Jv_31 = a4*sin(k)*sin(q1)*sin(q2)*sin(q4) - d7*sig_13 - a4*cos(q4)*sig_62 - a2*cos(q2)*sin(k)*sin(q1) - a3*cos(q1)*sin(k)*sin(q3) - d3*sin(k)*sin(q1)*sin(q2) - a3*cos(q2)*cos(q3)*sin(k)*sin(q1) - d5*sig_30;

    Jv_12 = d3*sig_66 + d5*sig_28 + d7*sig_12 + a2*cos(q2)*sin(k) + a3*cos(q3)*sig_61 - a4*sin(q4)*sig_66 - a2*cos(k)*cos(q1)*sin(q2) + a4*cos(q3)*cos(q4)*sig_61;
    Jv_22 = -sin(q1)*(a2*sin(q2) - d3*cos(q2) - d5*cos(q2)*cos(q4) + a3*cos(q3)*sin(q2) + a4*cos(q2)*sin(q4) - d7*cos(q2)*cos(q4)*cos(q6) + a4*cos(q3)*cos(q4)*sin(q2) + d5*cos(q3)*sin(q2)*sin(q4) + d7*cos(q3)*cos(q6)*sin(q2)*sin(q4) + d7*cos(q2)*cos(q5)*sin(q4)*sin(q6) - d7*sin(q2)*sin(q3)*sin(q5)*sin(q6) + d7*cos(q3)*cos(q4)*cos(q5)*sin(q2)*sin(q6));
    Jv_32 = a4*sin(q4)*sig_65 - d5*sig_26 - d7*sig_11 - a2*cos(k)*cos(q2) - a3*cos(q3)*sig_60 - d3*sig_65 - a2*cos(q1)*sin(k)*sin(q2) - a4*cos(q3)*cos(q4)*sig_60;

    Jv_13 = -d7*sig_18 - a4*cos(q4)*sig_55 - d5*sin(q4)*sig_55 - a3*sin(q3)*sig_66 - a3*cos(k)*cos(q3)*sin(q1);
    Jv_23 = d5*sin(q4)*sig_54 - d7*sig_16 + a3*cos(q1)*cos(q3) + a4*cos(q4)*sig_54 - a3*cos(q2)*sin(q1)*sin(q3);
    Jv_33 = d7*sig_19 + a4*cos(q4)*sig_56 + d5*sin(q4)*sig_56 + a3*sin(q3)*sig_65 - a3*cos(q3)*sin(k)*sin(q1);

    Jv_14 = d7*sig_21 + d5*sig_44 - a4*sin(q4)*sig_58 + a4*cos(q4)*sig_61;
    Jv_24 = d7*sig_20 + d5*sig_45 - a4*sin(q4)*sig_59 - a4*cos(q4)*sin(q1)*sin(q2);
    Jv_34 = a4*sin(q4)*sig_57 - d5*sig_43 - d7*sig_22 - a4*cos(q4)*sig_60;

    Jv_15 = -d7*sin(q6)*sig_8;
    Jv_25 = -d7*sin(q6)*sig_10;
    Jv_35 = d7*sin(q6)*sig_7;

    Jv_16 = -d7*sig_3;
    Jv_26 = -d7*sig_6;
    Jv_36 = d7*sig_2;

    // Sección que representa la velocindad angular
    Jw_41 = sig_10*(cos(q5)*sig_49 - sin(q5)*sig_48) - sig_13*sig_9 + (cos(q6)*sig_31 - sin(q6)*sig_30)*sig_6;
    Jw_51 = -(cos(q5)*sig_51 - sin(q5)*sig_50)*sig_7 - sig_14*sig_4 - (cos(q6)*sig_33 - sin(q6)*sig_32)*sig_2;
    Jw_61 = sig_5*sig_15 - sig_3*(sin(q6)*sig_34 - cos(q6)*sig_35) - (sin(q5)*sig_52 - cos(q5)*sig_53)*sig_8;

    Jw_42 = sig_10*(sin(q5)*sig_46 - cos(q5)*sin(q3)*sig_60) - sig_9*sig_11 - sig_6*(sin(q6)*sig_26 + cos(q6)*sig_27);
    Jw_52 = (sin(q5)*sig_47 - cos(q5)*sin(q3)*sig_61)*sig_7 - (sin(q6)*sig_28 + cos(q6)*sig_29)*sig_2 + sig_12*sig_4;
    Jw_62 = sig_3*(sin(q6)*sig_1 + cos(q6)*sig_17) - sig_5*(cos(q6)*sig_1 - sin(q6)*sig_17) - (sin(q5)*sig_37 + cos(q5)*sin(q1)*sin(q2)*sin(q3))*sig_8;

    Jw_43 = sig_9*sig_19 - sig_10*(cos(q5)*sig_57 - cos(q4)*sin(q5)*sig_56) - sig_6*(cos(q6)*sig_39 - sin(q4)*sin(q6)*sig_56);
    Jw_53 = -sig_18*sig_4 - (cos(q6)*sig_38 - sin(q4)*sin(q6)*sig_55)*sig_2 - (cos(q5)*sig_58 - cos(q4)*sin(q5)*sig_55)*sig_7;
    Jw_63 = (cos(q5)*sig_59 + cos(q4)*sin(q5)*sig_54)*sig_8 + (cos(q6)*sig_36 + sin(q4)*sin(q6)*sig_54)*sig_3 + sig_16*sig_5;

    Jw_44 = sin(q5)*sig_42*sig_10 - sig_6*(sin(q6)*sig_43 + cos(q5)*cos(q6)*sig_42) - sig_9*sig_22;
    Jw_54 = sig_21*sig_4 - (sin(q6)*sig_44 + cos(q5)*cos(q6)*sig_41)*sig_2 + sin(q5)*sig_41*sig_7;
    Jw_64 = (sin(q6)*sig_45 + cos(q5)*cos(q6)*sig_40)*sig_3 - sig_20*sig_5 - sin(q5)*sig_40*sig_8;

    Jw_45 = sig_10*sig_23 - cos(q6)*sig_6*sig_7 + sin(q6)*sig_9*sig_7;
    Jw_55 = sig_24*sig_7 - cos(q6)*sig_8*sig_2 - sin(q6)*sig_8*sig_4;
    Jw_65 = sig_25*sig_8 + cos(q6)*sig_10*sig_3 + sin(q6)*sig_10*sig_5;

    Jw_46 = sig_9*sig_2 + sig_6*sig_4;
    Jw_56 = sig_5*sig_2 - sig_3*sig_4;
    Jw_66 = sig_9*sig_3 + sig_6*sig_5;

    double J[6][6] = {{Jv_11,Jv_12,Jv_13,Jv_14,Jv_15,Jv_16},
                      {Jv_21,Jv_22,Jv_23,Jv_24,Jv_25,Jv_26},
                      {Jv_31,Jv_32,Jv_33,Jv_34,Jv_35,Jv_36},
                      {Jw_41,Jw_42,Jw_43,Jw_44,Jw_45,Jw_46},
                      {Jw_51,Jw_52,Jw_53,Jw_54,Jw_55,Jw_56},
                      {Jw_61,Jw_62,Jw_63,Jw_64,Jw_65,Jw_66}};
        
    cout << "Jacobiano :" << endl;
    for(int x=0;x<6;x++)  
    {
      for(int y=0;y<6;y++) 
        {
          cout << setw(10) << setprecision(4) << fixed << J[x][y]; // setw() = dimensión de la salida, fixed() y setprecision() = decimales
        }
    cout<<endl;
    }
}

int main ()
{
    jacob_geo(q1,q2,q3,q4,q5,q6,a2,d3,a3,a4,d5,d7,k);
    return 0;
}