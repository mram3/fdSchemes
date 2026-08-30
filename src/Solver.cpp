#include "Solver.h"

void Solver::computeFirstDeriv(
    Grid& mesh, Accuracy order
)
{
    if(order == Accuracy::firstOrder){
        int offset = 1; //for 1&2 order of accuracy

        mesh.f_phi.assign(mesh.N+1, 0.0);

        FDFA fd; BDFA bd;
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

        FDSA fd; CDSA cd; BDSA bd;
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

void Solver::computeSecondDeriv(
    Grid& mesh, Accuracy order
)
{
    if(order == Accuracy::firstOrder){
        int left_off = 1;
        int right_off = 1;
        FDFA fd; CDFA cd; BDFA bd;

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
        FDSA fd; CDSA cd; BDSA bd;

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

void Solver::implicitscheme(
    Grid& mesh, Accuracy order
)
{
    impSchemes pade;

    if(order == Accuracy::fourthOrder){
        mesh.f_phi.assign(mesh.N+1, 0.0);
        pade.pade4(mesh, mesh.f_phi);
    }
}