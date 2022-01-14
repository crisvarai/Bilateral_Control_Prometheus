#include <iostream>
#include <cmath>     // To use "pi" and "trigonometric" functions 
#include <iomanip>   // To "setw" and "setprecision" functions

using namespace std;

//----------------------- Declaración de variables del Controlador --------------------------//

// Variables para estimación de velocidad
double Ts; // Tiempo de muestreo
double xl, yl, zl, nl, bl1, bl2, bl3; // Variables de pose para humano
double xr, yr, zr, nr, br1, br2, br3; // Variables de pose para robot
double g1 = 1, g2 = 1, g3 = 1, g4 = 1, g5 = 1, g6 = 1; // Terminos de compensacion de gravedad
double yxpr, yypr, yzpr, ynpr, yb1pr, yb2pr, yb3pr; // Dinamica del filtro
static double yxr = 0, yyr = 0, yzr = 0, ynr = 0, yb1r = 0, yb2r = 0, yb3r = 0; // Integral del filtro
double kyi = 1; // Ganancia del FILTRO

// Variables de fuerza
double f1, f2, f3, f4, f5, f6; // Vector de fuerza definido en Ec.(15)

// Elementos del Jacobiano
double J11, J12, J13, J14, J15, J16;
double J21, J22, J23, J24, J25, J26;
double J31, J32, J33, J34, J35, J36;
double J41, J42, J43, J44, J45, J46;
double J51, J52, J53, J54, J55, J56;
double J61, J62, J63, J64, J65, J66;

// Definimos las variables de torque
double tau1, tau2, tau3, tau4, tau5, tau6;

// Definimos la corriente a cada motor
double id1, id2, id3, id4, id5, id6;

// Definimos corriente Maxima
double IM = 0.5;

// Ganancias del control
double kr1 = 1, kr2 = 1, kr3 = 1, kr4 = 1, kr5 = 1, kr6 = 1; // Ganancia Proporcional
double kd1 = 1, kd2 = 1, kd3 = 1, kd4 = 1, kd5 = 1, kd6 = 1; // Ganancia Derivativa

// Variables no definidas, obtención de velocidad (Preguntar a Guerrero)
double x, y, z, ybr1, ybr2, ybr3;

// Variables no definidas, vector de fuerza (Preguntar a Guerrero)
double dr1, dr2, dr3, dr4, dr5, dr6;

int main() {
    
    // Código para obtener la velocidad de las variables del sensor
    yxpr = x-kyi-yxr;
    yxr += Ts*yxpr;

    yypr = y-kyi*yyr;
    yyr += Ts*yypr;

    yzpr = z-kyi*yzr;
    yzr += Ts*yzpr;

    ynpr = nr-kyi*ynr;
    ynr += Ts*ynpr;

    yb1pr = br1-kyi*yb1r; // revisar a detalle
    ybr1 += Ts*yb1pr;

    yb2pr = br2-kyi*yb2r; // revisar a detalle
    ybr2 += Ts*yb2pr;

    yb3pr = br3-kyi*yb3r; // revisar a detalle
    ybr3 += Ts*yb3pr;

    // Definimos cada elemento del jacobiano que está en función del dato
    // "q" dado por el encoder
    // J11 = cos(q1)+sin(q1); \\ Por ejemplo
    // ...
    // J66 = función de "q";

    // Definición de los elementos del vector de fuerza
    f1 = -kr1*(xr-xl)-dr1+yxpr;
    f2 = -kr2*(yr-yl)-dr2*yypr;
    f3 = -kr3*(zr-zl)-dr3*yzpr;
    f4 = -kr4*(nl*br1-nr*bl1-bl3*br2-bl2*br3)-0.5*dr4*(br2*yb3pr-br3*yb2pr-br1*ynpr)+nr*yb1pr; // Revisar con Guerrero, problema en los "()"
    f5 = -kr5*(nl*br2-nr*bl2-bl1*br3-bl3*br1)-0.5*dr5*(br1*yb3pr-br3*yb1pr+br2*ynpr)-nr*yb2pr;
    f6 = -kr6*(nl*br3-nr*bl3-bl2*br1-bl1*br2)-0.5*dr6*(br1*yb2pr-br2*yb1pr-br3*ynpr)+nr*yb3pr;

    // Calculamos el torque para cada actuador
    tau1 = J11*f1+J21*f2+J31*f3+J41*f4+J51*f5+J61*f6+g1;
    tau2 = J12*f1+J22*f2+J32*f3+J42*f4+J52*f5+J62*f6+g2;
    tau3 = J13*f1+J23*f2+J33*f3+J43*f4+J53*f5+J63*f6+g3;
    tau4 = J14*f1+J24*f2+J34*f3+J44*f4+J54*f5+J64*f6+g4;
    tau5 = J15*f1+J25*f2+J35*f3+J45*f4+J55*f5+J65*f6+g5;
    tau6 = J16*f1+J26*f2+J36*f3+J46*f4+J56*f5+J66*f6+g6;

    // Convertimos el torque en corriente
    // id1 = tau1/km;
    // id2 = tau2/km;
    // id3 = tau3/km;
    // id4 = tau4/km;
    // id5 = tau5/km;
    // id6 = tau6/km;

    // Saturamos por corriente para proteger el motor
    if(id1 > IM) { 
        id1 = IM; 
    }
    if(id1 < -IM) { 
        id1 = -IM; 
    }
    if(id2 > IM) { 
        id2 = IM; 
    }
    if(id2 < -IM) { 
        id2 = -IM; 
    }
    if(id3 > IM) { 
        id3 = IM; 
    }
    if(id3 < -IM) { 
        id3 = -IM; 
    }
    if(id4 > IM) { 
        id4 = IM; 
    }
    if(id4 < -IM) { 
        id4 = -IM; 
    }
    if(id5 > IM) { 
        id5 = IM; 
    }
    if(id5 < -IM) { 
        id5 = -IM; 
    }
    if(id6 > IM) { 
        id6 = IM; 
    }
    if(id6 < -IM) { 
        id6 = -IM; 
    }

    // cout << setprecision(6) << fixed << "yxpr :" << yxpr << endl;

    return 0;
}


