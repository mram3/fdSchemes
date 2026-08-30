#include "ImplicitSchemes.h"
#include "MathTools.h"

#include <vector>

using namespace std;

void impSchemes::pade4(
    const Grid& mesh, vector<double>& results
)const{
    int N = mesh.N;
    double dx = mesh.x[1] - mesh.x[0];
    vector<double> D(N+1);
    vector<double> E(N+1);
    vector<double> C(N+1);
    vector<double> b(N+1);

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