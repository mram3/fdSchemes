#ifndef EXPLICITSCHEMES_H
#define EXPLICITSCHEMES_H

#include "Grid.h"

class expSchemes //PARENT
{
public:

    virtual ~expSchemes() = default;//virtual destructor

    virtual double firstDeriv(
        const Grid& mesh, int i
    )const = 0; // = 0 -> pure virtual function

    virtual double secondDeriv(
        const Grid& mesh, int i
    )const = 0;
};

/*Abbreviations: 
        CD - Centre Difference   1A - First order accurate
        FD - Forward Difference  2A - Second order accurate
        BD - Backward Difference 4A - Fourth order accurate
*/

class CD4A: public expSchemes
{
public:
    double firstDeriv(
       const Grid& mesh, int i
    ) const override;

    double secondDeriv(
        const Grid& mesh, int i
    ) const override;
};
class CD2A: public expSchemes //CHILD 1
{
public:
    double firstDeriv(
        const Grid& mesh, int i
    ) const override;

    double secondDeriv(
       const Grid& mesh, int i
    ) const override;
};

class CD1A: public expSchemes //CHILD 2
{
public:
    double firstDeriv(
        const Grid& mesh, int i
    ) const override;

    double secondDeriv(
        const Grid& mesh, int i
    ) const override;
};

class FD4A: public expSchemes
{
public:
    double firstDeriv(
       const Grid& mesh, int i
    ) const override;

    double secondDeriv(
        const Grid& mesh, int i
    ) const override;
};
class FD2A: public expSchemes //CHILD 3
{
public:
    double firstDeriv(
        const Grid& mesh, int i
    ) const override;

    double secondDeriv(
        const Grid& mesh, int i
    ) const override;
};

class FD1A: public expSchemes //CHILD 4
{
public:
    double firstDeriv(
        const Grid& mesh, int i
    ) const override;

    double secondDeriv(
        const Grid& mesh, int i
    ) const override;
};

class BD4A: public expSchemes
{
public:
    double firstDeriv(
       const Grid& mesh, int i
    ) const override;

    double secondDeriv(
        const Grid& mesh, int i
    ) const override;
};
class BD2A: public expSchemes //CHILD 5
{
public:
    double firstDeriv(
        const Grid& mesh, int i
    ) const override;

    double secondDeriv(
        const Grid& mesh, int i
    ) const override;
};

class BD1A: public expSchemes //CHILD 6
{
public:
    double firstDeriv(
       const Grid& mesh, int i
    ) const override;

    double secondDeriv(
        const Grid& mesh, int i
    ) const override;
};

#endif