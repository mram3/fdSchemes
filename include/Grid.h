#ifndef GRID_H
#define GRID_H

#include "Type.h"
#include <vector>

class Grid
{
public: 

    int N; //Number of grid points
    double xmin, xmax; //dimensions
    std::vector <double> x, phi, f_phi, s_phi;

    Grid();

    Grid(
        int N,
        double xmin,
        double xmax,
        MeshType type
    );

    void Phi();
};

#endif