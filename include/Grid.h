#ifndef GRID_H
#define GRID_H

#include "Type.h"
#include <vector>
#include <functional>

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
        std::function<double(int, double, double, int)> meshFunction
    );

    void Phi(std::function<double(double)> exactFunction);
};

#endif