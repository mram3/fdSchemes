#ifndef MATHTOOLS_H
#define MATHTOOLS_H

#include<vector>

class MathTools{
public:
    static void TDMA(
        const std::vector<double>& C,
        const std::vector<double>& D,
        const std::vector<double>& E,
        const std::vector<double>& b,
        std::vector<double>& results
    );
};
#endif