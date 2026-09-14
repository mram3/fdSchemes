#ifndef SOLVER_H
#define SOLVER_H

#include "Grid.h"
#include "ExplicitSchemes.h"
#include "Type.h"
#include "ImplicitSchemes.h"
#include "MathTools.h"

class Solver{
public:

    void explicitFirstDeriv(
        Grid& mesh, Accuracy order
    );

    void explicitSecondDeriv(
        Grid& mesh, Accuracy order
    );

    void implicitFirstDeriv(
        Grid& mesh, Accuracy order
    );

    void implicitSecondDeriv(
        Grid& mesh, Accuracy order
    );

};

#endif