/*
    Program to get desired derivative of desired order

Compiling Instruction: g++ 1_PolynomialDS.cpp src/*.cpp -Iinclude -std=c++17 -o pds && ./pds
*/

#include <iostream>
#include <cmath>
#include <fstream>
#include <string>

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

double exactSecondDeriv(double x){
    return -2*exp(-x)*cos(x);
}

int main(){
    int N, A, O;
    double xmax, xmin;
    string str;

    cout << "=================================================\n";
    cout << "        POLYNOMIAL DIFFERENCING SCHEMES\n";
    cout << "=================================================\n";

    cout << "\nEnter the number of node points: ";
    cin >> N;
    cout << endl;

    cout << "Enter xmin: ";
    cin >> xmin;
    cout << endl;

    cout << "Enter xmax: ";
    cin >> xmax;
    cout << endl;

    Grid mesh(N, xmin, xmax, meshFunction);
    mesh.Phi(exactFunction);
    cout << str<<"Mesh is created successfully\n\n";

    cout << "Desired order of derivative\n";
    cout << "Enter 1 for first derivative\n";
    cout << "Enter 2 for second derivative\n";
    cout << "Enter : ";
    cin >> O;
    cout << endl;

    if(O != 1 && O != 2){
        cout << "\nError: Enter only 1st or 2nd order derivative\n";
        return 1;
    }

    cout << "Desired order of accuracy\n";
    cout << "Enter 1 for first order\n";
    cout << "Enter 2 for second order\n";
    cout << "Enter: ";
    cin >> A;
    cout << endl;

    Accuracy accuracy;
    if(A==1){
        accuracy = Accuracy::firstOrder;
        str = "First order accurate ";
    }
    else if(A==2){
        accuracy = Accuracy::secondOrder;
        str = "Second order accurate ";
    }
    else{
        cout<<"\nError: Enter only 1 or 2. Program Terminated\n";
        return 1;
    }

    //Storing the outputs in a .csv file
    ofstream file ("1_PolynomialDS_Output.csv");

    Solver solve;

    if(O==1){
        solve.explicitFirstDeriv(mesh, accuracy);
        //Writing the output file to csv file
        file <<str<<"First Derivative"<<endl;
        file<<"x,numerical solution, exact solution, error"<<endl;
        for(int i = 0; i <= N; i++){
            file << mesh.x[i] << "," << mesh.f_phi[i] << "," 
            << exactFirstDeriv(mesh.x[i]) << "," 
            << fabs(mesh.f_phi[i]-exactFirstDeriv(mesh.x[i])) << endl; 
        }
    }
    else if(O==2){
        solve.explicitSecondDeriv(mesh, accuracy);
        //Writing the output file to csv file
        file<<str<<"Second Derivative"<<endl;
        file<<"x,numerical solution, exact solution, error"<<endl;
        for(int i = 0; i <= N; i++){
            file << mesh.x[i] << "," << mesh.s_phi[i] << "," 
            << exactSecondDeriv(mesh.x[i]) << "," 
            << fabs(mesh.s_phi[i]-exactSecondDeriv(mesh.x[i])) << endl; 
        }
    }

    cout << "Derivative values are stored in 1_PolynomialDS_Output.csv"<<endl;
    return 0;
}