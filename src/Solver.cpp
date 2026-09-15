#include "Solver.h"

void Solver::explicitFirstDeriv(
    Grid& mesh, Accuracy order
)
{
    if(order == Accuracy::firstOrder){
        int offset = 1; //for 1&2 order of accuracy

        mesh.f_phi.assign(mesh.N+1, 0.0);

        FD1A fd; BD1A bd;
        for(int i = 0; i < offset; i++){
            mesh.f_phi[i] = fd.firstDeriv(mesh, i);
        }

        for(int i = offset; i<=mesh.N; i++){
            mesh.f_phi[i] = bd.firstDeriv(mesh, i);
        }
    }

    else if(order == Accuracy::secondOrder){
        int offset = 1;

        mesh.f_phi.assign(mesh.N+1, 0.0);

        FD2A fd; CD2A cd; BD2A bd;
        for(int i = 0; i < offset; i++){
            mesh.f_phi[i] = fd.firstDeriv(mesh, i);
        }

        for(int i = offset; i<mesh.N-offset; i++){
            mesh.f_phi[i] = cd.firstDeriv(mesh, i);
        }

        for(int i = mesh.N-offset; i<=mesh.N; i++){
            mesh.f_phi[i] = bd.firstDeriv(mesh, i);
        }
    }
}

void Solver::explicitSecondDeriv(
    Grid& mesh, Accuracy order
)
{
    if(order == Accuracy::firstOrder){
        int left_off = 1;
        int right_off = 1;
        FD1A fd; CD1A cd; BD1A bd;

        mesh.s_phi.assign(mesh.N+1, 0.0);

        for(int i = 0; i < left_off; i++){
            mesh.s_phi[i] = fd.secondDeriv(mesh, i);
        }

        for(int i = left_off; i<mesh.N-right_off; i++){
            mesh.s_phi[i] = cd.secondDeriv(mesh, i);
        }
        for(int i = mesh.N-right_off; i<= mesh.N; i++){
            mesh.s_phi[i] = bd.secondDeriv(mesh, i);
        }
    }

    else if(order == Accuracy::secondOrder){
        int left_off = 1;
        int right_off = 2;
        FD2A fd; CD2A cd; BD2A bd;

        mesh.s_phi.assign(mesh.N+1, 0.0);

        for(int i = 0; i < left_off; i++){
            mesh.s_phi[i] = fd.secondDeriv(mesh, i);
        }

        for(int i = left_off; i<mesh.N-right_off; i++){
            mesh.s_phi[i] = cd.secondDeriv(mesh, i);
        }
        for(int i = mesh.N-right_off; i<= mesh.N; i++){
            mesh.s_phi[i] = bd.secondDeriv(mesh, i);
        }
    }
}

void Solver::implicitFirstDeriv(
    Grid& mesh, Accuracy order
)
{
    if(order == Accuracy::fourthOrder){
        pade4 pade;
        mesh.f_phi.assign(mesh.N+1, 0.0);
        pade.firstDeriv(mesh, mesh.f_phi);
    }

    else if(order == Accuracy::sixthOrder){
        pade6 pade;
        mesh.f_phi.assign(mesh.N+1, 0.0);
        pade.firstDeriv(mesh, mesh.f_phi);
    }
}