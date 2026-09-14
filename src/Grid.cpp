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
    function<double(int, double, double, int)> meshFunction
)
{
    N = N_;
    xmin = xmin_;
    xmax = xmax_;
    x.assign(N+1, 0.0);

    for(int i = 0; i <= N; i++){
        x[i] = meshFunction(i, xmin, xmax, N);
    }
}

void Grid::Phi(function<double(double)> exactFunction)
{
    phi.assign(N+1, 0.0);

    for(int i = 0; i <= N; i++){
        phi[i] = 
            exactFunction(x[i]);
    }
}