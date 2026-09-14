/*
Compiling Instruction: g++ 2_CCDS.cpp src/*.cpp -Iinclude -std=c++17 -o ccds && ./ccds
*/
#include <iostream>
#include <cmath>
#include "vector"
#include <fstream>

#include "Solver.h"


using namespace std;

double meshFunction(int i, double xmin, double xmax, int N){
    double PI = acos(-1);
    return -((xmax-xmin)/2)*cos(i*PI/N) + ((xmax+xmin)/2);
}

double exactFunction(double x){
    return exp(-x)*sin(x);
}

double exactFirstDeriv(double x){
    return exp(-x) * (cos(x)-sin(x));
}

int main(){

    int N = 100;
    double xmin = 0;
    double xmax = 10;
    vector<double> A, B, C;
    Grid mesh(N, xmin, xmax, meshFunction);
    mesh.Phi(exactFunction);

    Solver solve;
    solve.implicitFirstDeriv(mesh, Accuracy::fourthOrder);
    A = mesh.f_phi;

    solve.explicitFirstDeriv(mesh, Accuracy::secondOrder);
    B = mesh.f_phi;
    for(int i = 0; i<=N;i++){
        C.push_back(exactFirstDeriv(mesh.x[i]));
    }

    ofstream file("Compact_vs_ConventionalFD.csv");

    file << "Performance of Compact Vs Conventional FD using same stencil sizes (3)"<<endl;
    file << "x," << "Pade4,"<<"SecondOrder FD," << "Exact"<<endl;
    for(int i = 0; i<=N;i++){
        file << mesh.x[i] << "," << A[i] << "," << B[i] << "," << C[i]<<endl;
    }

    cout << "Results are stored in a Compact_vs_ConventionalFD.csv file";
    return 0;
}