#ifndef IMPLICITSCHEMES_H
#define IMPLICITSCHEMES_H

#include "Grid.h"

class impSchemes //PARENT
{
public:

    virtual ~impSchemes() = default;//virtual destructor

    virtual void firstDeriv(
        const Grid& mesh, std::vector<double>& results
    )const = 0; // = 0 -> pure virtual function

    virtual void secondDeriv(
        const Grid& mesh, std::vector<double>& results
    )const = 0;
};

//Central Compact Finite Difference schemes
class pade4 : public impSchemes
{
public:
    void firstDeriv(
        const Grid& mesh, std::vector<double>& results
    ) const override;

    void secondDeriv(
        const Grid& mesh, std::vector<double>& results
    ) const override;
};

class pade6 : public impSchemes
{
public:
    void firstDeriv(
        const Grid& mesh, std::vector<double>& results
    ) const override;

    void secondDeriv(
        const Grid& mesh, std::vector<double>& results
    ) const override;
};
/*class impSchemes{
public:
    void pade4(
        const Grid& mesh, std::vector<double>& results
    )const; //cental compact finite difference

    void ucfd(Grid& mesh) const; //upwind compact finite difference
};*/

#endif