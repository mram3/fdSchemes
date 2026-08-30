#ifndef EXPLICITSCHEMES_H
#define EXPLICITSCHEMES_H

#include "Grid.h"

class Schemes //PARENT
{
public:

    virtual ~Schemes() = default;//virtual destructor

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

class CDSA: public Schemes //CHILD 1
{
public:
    double firstDeriv(
        const Grid& mesh, int i
    ) const override;

    double secondDeriv(
       const Grid& mesh, int i
    ) const override;
};

class CDFA: public Schemes //CHILD 2
{
public:
    double firstDeriv(
        const Grid& mesh, int i
    ) const override;

    double secondDeriv(
        const Grid& mesh, int i
    ) const override;
};

class FDSA: public Schemes //CHILD 3
{
public:
    double firstDeriv(
        const Grid& mesh, int i
    ) const override;

    double secondDeriv(
        const Grid& mesh, int i
    ) const override;
};

class FDFA: public Schemes //CHILD 4
{
public:
    double firstDeriv(
        const Grid& mesh, int i
    ) const override;

    double secondDeriv(
        const Grid& mesh, int i
    ) const override;
};

class BDSA: public Schemes //CHILD 5
{
public:
    double firstDeriv(
        const Grid& mesh, int i
    ) const override;

    double secondDeriv(
        const Grid& mesh, int i
    ) const override;
};

class BDFA: public Schemes //CHILD 6
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