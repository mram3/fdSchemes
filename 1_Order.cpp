/*
Compiling Instruction: g++ 1_Order.cpp src/*.cpp -Iinclude -std=c++17 -o order && ./order
*/

#include "Solver.h"
#include <cmath>
#include <fstream>
#include <vector>
#include <iostream>

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

double exactSecondDeriv(double x){
    return -2*exp(-x)*cos(x);
}

double getError_n_write(
    int N, double xmin, double xmax, int O, int A, ofstream& file
)
{

    Grid mesh (N, xmin, xmax, meshFunction);
    mesh.Phi(exactFunction);
    Solver solve;
    double err, maxErr = 0;
    if(O==1){
        if(A==1){
            solve.explicitFirstDeriv(mesh, Accuracy::firstOrder);
            file << "FIRST ORDER ACCURATE FIRST DERIVATIVE" << endl;
            file << "N = " << N << endl;
            file << "x," << "exact," << "numerical," << "error" << endl;
            for(int i = 0; i<=N; i++){
                err = abs(mesh.f_phi[i]-exactFirstDeriv(mesh.x[i]));
                if(err > maxErr) maxErr =  err;
                file<< mesh.x[i] << ","
                    << exactFirstDeriv(mesh.x[i]) << ","
                    << mesh.f_phi[i] << ","
                    << err << endl;
            }
            file << endl;
            return maxErr;
        }
        else if(A==2){
            solve.explicitFirstDeriv(mesh, Accuracy::secondOrder);
            file << "SECOND ORDER ACCURATE FIRST DERIVATIVE" << endl;
            file << "N = " << N << endl;
            file << "x," << "exact," << "numerical," << "error" << endl;
            for(int i = 0; i<=N; i++){
                err = abs(mesh.f_phi[i]-exactFirstDeriv(mesh.x[i]));
                if(err > maxErr) maxErr =  err;
                file<< mesh.x[i] << ","
                    << exactFirstDeriv(mesh.x[i]) << ","
                    << mesh.f_phi[i] << ","
                    << err << endl;
            }
            file << endl;
            return maxErr;
        }
    }
    else if(O==2){
        if(A==1){
            solve.explicitSecondDeriv(mesh, Accuracy::firstOrder);
            file << "FIRST ORDER ACCURATE SECOND DERIVATIVE" << endl;
            file << "N = " << N << endl;
            file << "x," << "exact," << "numerical," << "error" << endl;
            for(int i = 0; i<=N; i++){
                err = abs(mesh.s_phi[i]-exactSecondDeriv(mesh.x[i]));
                if(err > maxErr) maxErr =  err;
                file<< mesh.x[i] << ","
                    << exactSecondDeriv(mesh.x[i]) << ","
                    << mesh.s_phi[i] << ","
                    << err << endl;
            }
            file << endl;
            return maxErr;
        }
        else if(A==2){
            solve.explicitSecondDeriv(mesh, Accuracy::secondOrder);
            file << "SECOND ORDER ACCURATE SECOND DERIVATIVE" << endl;
            file << "N = " << N << endl;
            file << "x," << "exact," << "numerical," << "error" << endl;
            for(int i = 0; i<=N; i++){
                err = abs(mesh.s_phi[i]-exactSecondDeriv(mesh.x[i]));
                if(err > maxErr) maxErr =  err;
                file<< mesh.x[i] << ","
                    << exactSecondDeriv(mesh.x[i]) << ","
                    << mesh.s_phi[i] << ","
                    << err << endl;
            }
            file << endl;
            return maxErr;
        }
    }

    return 0.0;
}

int main(){


    int N1 = 50;
    int N2 = 100;

    double xmin = 0;
    double xmax = 10;

    //First Derivative
    ofstream file1("First_Derivative.csv");

    cout << "========================================" << endl;
    cout << "           FIRST DERIVATIVE \n";
    cout << "========================================" << endl;
    
    //First Order
    double E1 = getError_n_write(N1, xmin,xmax,1,1, file1);
    double E2 = getError_n_write(N2, xmin,xmax,1,1, file1);

    double accuracy1 = log(E1/E2)/log((double)N2/N1);

    cout << "                First Order\n";
    cout << "Maximum Error when N = " << N1 << " is " << E1 << endl;
    cout << "Maximum Error when N = " << N2 << " is " << E2 << endl;
    cout << " Observed Order of Accuracy = " << accuracy1 << endl;

    //Second Order
    E1 = getError_n_write(N1, xmin,xmax,1,2, file1);
    E2 = getError_n_write(N2, xmin,xmax,1,2, file1);

    double accuracy2 = log(E1/E2)/log((double)N2/N1);

    cout << "                Second Order\n";
    cout << "Maximum Error when N = " << N1 << " is " << E1 << endl;
    cout << "Maximum Error when N = " << N2 << " is " << E2 << endl;
    cout << " Observed Order of Accuracy = " << accuracy2 << endl;

    //Second Derivative 

    ofstream file2("Second_Derivative.csv");

    cout << "========================================" << endl;
    cout << "           SECOND DERIVATIVE \n";
    cout << "========================================" << endl;
    
    //First Order
    E1 = getError_n_write(N1, xmin,xmax,2,1, file2);
    E2 = getError_n_write(N2, xmin,xmax,2,1, file2);

    double accuracy3 = log(E1/E2)/log((double)N2/N1);

    cout << "                First Order\n";
    cout << "Maximum Error when N = " << N1 << " is " << E1 << endl;
    cout << "Maximum Error when N = " << N2 << " is " << E2 << endl;
    cout << " Observed Order of Accuracy = " << accuracy3 << endl;

    //Second Order
    E1 = getError_n_write(N1, xmin,xmax,2,2, file2);
    E2 = getError_n_write(N2, xmin,xmax,2,2, file2);

    double accuracy4 = log(E1/E2)/log((double)N2/N1);

    cout << "                Second Order\n";
    cout << "Maximum Error when N = " << N1 << " is " << E1 << endl;
    cout << "Maximum Error when N = " << N2 << " is " << E2 << endl;
    cout << " Observed Order of Accuracy = " << accuracy4 << endl;

    cout << "========================================" << endl;
    cout << "All data stored in First_Derivative.csv & Second_Derivative.csv files";
    return 0;
}