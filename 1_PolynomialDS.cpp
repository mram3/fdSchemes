/*
    Program to get desired derivative of desired order
*/
#include <iostream>
#include "Grid.h"
#include "ExplicitSchemes.h"
#include "Solver.h"
#include "Type.h"

using namespace std;

int main(){
    int N, A, O;
    double xmax, xmin;

    cout << "=================================================\n";
    cout << "        POLYNOMIAL DIFFERENCING SCHEMES\n";
    cout << "=================================================\n";

    cout << "Enter the number of node points: ";
    cin >> N;
    cout << endl;

    cout << "Enter xmax: ";
    cin >> xmax;
    cout << endl;

    cout << "Enter xmax: ";
    cin >> xmax;
    cout << endl;

    cout << "Enter desired order of derivatives" << endl;
    cout << "Enter 1 for first derivative" << endl;
    cout << "Enter 2 for second derivative" << endl;
    cout << "Enter : ";
    cin >> O;
    cout << endl;

    cout << "Enter desired order of accuracy: ";
    cin >> A;
    cout << endl;

    Grid mesh( //initialises the grid
        N, xmin, xmax, MeshType::nonUniform
    );

    mesh.Phi(); //sets phi values 

    Solver solve;
    
}