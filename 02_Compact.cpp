/*

Program to implement the compact scheme and compare it with conventional FD scheme

1) meshFunction : Control the spacing between grids
2) exactFunction : Function for which derivatives are to be calculated
3) Results are stored in 02_Compact_results.csv
4) 02_Compact_plot.py : Plot the derivatives and errors to visually compare

Compiling Instruction: g++ 2_CCDS.cpp src/*.cpp -Iinclude -std=c++17 -o ccds && ./ccds
*/

#include <iostream>
#include <cmath>
#include "vector"
#include <fstream>

#include "Solver.h"

using namespace std;

double meshFunction(int i, double xmin, double xmax, int N){
    return xmin + i * (xmax-xmin)/(double)N; //uniform spacing
}

//Viscous Burger Equation (Shock profile)
// u * du/dx = nu * d^2(u)/dx^2
const double U_inf = 1.0;
const double a = 2.0; // a = U / (2*nu); nu = 0.25 (say)

double exactFunction(double x){
    return -U_inf*tanh(a*x);
}

//1st Deriv: u'(x) = -U*a * sech^2(a*x)
double exactFirstDeriv(double x){
    double cosh_ax = cosh(a * x);
    double sech2_ax = 1.0 / (cosh_ax * cosh_ax);
    
    return -U_inf * a * sech2_ax;
}

//2nd Deriv: u''(x) = 2U*a^2 * tanh(a*x) * sech^2(a*x)
double exactSecondDeriv(double x){
    double cosh_ax = cosh(a * x);
    double sech2_ax = 1.0 / (cosh_ax * cosh_ax);
    double tanh_ax = tanh(a * x);
    
    return 2 * U_inf * a * a * tanh_ax * sech2_ax;
}

int main(){

    int N = 100;
    double xmin = -5.0;
    double xmax = 5.0;
    vector<double> A, B, C;
    Grid mesh(N, xmin, xmax, meshFunction);
    mesh.Phi(exactFunction);

    ofstream file("02_Compact_results.csv");
    Solver solve;

    //First Derivative
    solve.implicitFirstDeriv(mesh, Accuracy::fourthOrder);
    A = mesh.f_phi;

    solve.explicitFirstDeriv(mesh, Accuracy::fourthOrder);
    B = mesh.f_phi;
    for(int i = 0; i<=N;i++){
        C.push_back(exactFirstDeriv(mesh.x[i]));
    }

    //write first derivative results
    file << "Performance of Compact Vs Conventional FD"<<endl;
    file << "First Derivative"<<endl;
    file << "x," << "Pade4,"<<"FourthOrder FD," << "Exact"<<endl;
    for(int i = 0; i<=N;i++){
        file << mesh.x[i] << "," << A[i] << "," << B[i] << "," << C[i]<<endl;
    }

    //second derivative
    A.clear(); B.clear(); C.clear();
    solve.implicitSecondDeriv(mesh, Accuracy::fourthOrder);
    A = mesh.s_phi;

    solve.explicitSecondDeriv(mesh, Accuracy::fourthOrder);
    B = mesh.s_phi;
    for(int i = 0; i<=N;i++){
        C.push_back(exactSecondDeriv(mesh.x[i]));
    }

    //write second derivative results
    file << "Second Derivative"<<endl;
    file << "x," << "Pade4,"<<"FourthOrder FD," << "Exact"<<endl;
    for(int i = 0; i<=N;i++){
        file << mesh.x[i] << "," << A[i] << "," << B[i] << "," << C[i]<<endl;
    }

    cout << "Results are stored in 02_Compact_results.csv file";
    return 0;
}