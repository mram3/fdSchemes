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
        CD - Centre Difference   FA - First order accurate
        FD - Forward Difference  SA - Second order accurate
        BD - Backward Difference
*/

class CDSA: public expSchemes //CHILD 1
{
public:
    double firstDeriv(
        const Grid& mesh, int i
    ) const override;

    double secondDeriv(
       const Grid& mesh, int i
    ) const override;
};

class CDFA: public expSchemes //CHILD 2
{
public:
    double firstDeriv(
        const Grid& mesh, int i
    ) const override;

    double secondDeriv(
        const Grid& mesh, int i
    ) const override;
};

class FDSA: public expSchemes //CHILD 3
{
public:
    double firstDeriv(
        const Grid& mesh, int i
    ) const override;

    double secondDeriv(
        const Grid& mesh, int i
    ) const override;
};

class FDFA: public expSchemes //CHILD 4
{
public:
    double firstDeriv(
        const Grid& mesh, int i
    ) const override;

    double secondDeriv(
        const Grid& mesh, int i
    ) const override;
};

class BDSA: public expSchemes //CHILD 5
{
public:
    double firstDeriv(
        const Grid& mesh, int i
    ) const override;

    double secondDeriv(
        const Grid& mesh, int i
    ) const override;
};

class BDFA: public expSchemes //CHILD 6
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