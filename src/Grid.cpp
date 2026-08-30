#include "Grid.h"

#include <vector>
#include <cmath>

using namespace std;

Grid::Grid(){
    N = 0;
    xmin = 0; xmax = 0;
}

Grid::Grid(
    int N_,
    double xmin_,
    double xmax_,
    MeshType type
)
{
    N = N_;
    xmin = xmin_;
    xmax = xmax_;
    x.assign(N+1, 0.0);

    if(type == MeshType::nonUniform){
        const double PI = acos(-1.0);
        for (int i = 0; i <= N; i++){
            x[i] =
                -((xmax-xmin)/2)*cos(i*PI/N) + ((xmax+xmin)/2);
        }
    }

    else if(type == MeshType::Uniform){
        double dy = (xmax-xmin)/N;
        for(int i = 0; i<=N; i++){
            x[i] = xmin + i*dy;
        }
    }
}

void Grid::Phi()
{
    phi.assign(N+1, 0.0);

    for(int i = 0; i <= N; i++){
        phi[i] = 
            exp(-x[i])*sin(x[i]);
    }
}