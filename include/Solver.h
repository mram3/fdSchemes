#ifndef SOLVER_H
#define SOLVER_H

#include "Grid.h"
#include "ExplicitSchemes.h"
#include "Type.h"
#include "ImplicitSchemes.h"
#include "MathTools.h"

class Solver{
public:

    void computeFirstDeriv(
        Grid& mesh, Accuracy order
    );

    void computeSecondDeriv(
        Grid& mesh, Accuracy order
    );

    void implicitscheme(
        Grid& mesh, Accuracy order
    );
};

#endif