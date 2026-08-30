#ifndef IMPLICITSCHEMES_H
#define IMPLICITSCHEMES_H

#include "Grid.h"

class impSchemes{
public:
    void pade4(
        const Grid& mesh, std::vector<double>& results
    )const; //cental compact finite difference

    void ucfd(Grid& mesh) const; //upwind compact finite difference
};

#endif