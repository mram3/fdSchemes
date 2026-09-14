/*
Implicit Schemes are written only for structured meshes
*/

#include "ImplicitSchemes.h"
#include "MathTools.h"

#include <vector>

using namespace std;

void pade4::firstDeriv(
    const Grid& mesh, vector<double>& results
)const 
{
    int N = mesh.N;
    double dx = mesh.x[1] - mesh.x[0];
    vector<double> D(N+1);//main diagonal 
    vector<double> E(N+1);//super diagonal (upper)
    vector<double> C(N+1);//sub diagonal (lower)
    vector<double> b(N+1);//RHS

    for(int i = 1; i <= N-1; i++){
        D[i] = 1;
        E[i] = 0.25;
        C[i] = 0.25;
        b[i] = 0.75*(mesh.phi[i+1]-mesh.phi[i-1])/dx;
    }
    D[0] = 1; D[N] = 1;
    E[0] = 2; E[N] = 0;
    C[0] = 0; C[N] = 2;
    b[0] = 0.5*(-5*mesh.phi[0]+4*mesh.phi[1]+mesh.phi[2])/dx;
    b[N] = 0.5*(5*mesh.phi[N]-4*mesh.phi[N-1]-mesh.phi[N-2])/dx;

    MathTools::TDMA(C, D, E, b, results);
}

void pade4::secondDeriv(
    const Grid& mesh, vector<double>& results
)const 
{
    int N = mesh.N;
    double dx = mesh.x[1] - mesh.x[0];
    vector<double> D(N+1);//main diagonal 
    vector<double> E(N+1);//super diagonal (upper)
    vector<double> C(N+1);//sub diagonal (lower)
    vector<double> b(N+1);//RHS

    for(int i = 1; i <= N-1; i++){
        D[i] = 1.0;
        E[i] = 0.1;
        C[i] = 0.1;
        b[i] = 6.0*(mesh.phi[i+1]+mesh.phi[i-1]
                  -2.0*mesh.phi[i])/(5.0*dx*dx);
    }
    D[0] = 1; D[N] = 1;
    E[0] = 11; E[N] = 0;
    C[0] = 0; C[N] = 11;
    b[0] = 13*mesh.phi[0]-27*mesh.phi[1]+15*mesh.phi[2]-mesh.phi[3];
    b[0] /= dx*dx;
    b[N] = 13*mesh.phi[N]-27*mesh.phi[N-1]+15*mesh.phi[N-2]-mesh.phi[N-3];
    b[N] /= dx*dx;

    MathTools::TDMA(C, D, E, b, results);
}