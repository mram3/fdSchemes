#include "MathTools.h"

using namespace std;

void MathTools::TDMA(
    const vector<double>& C,
    const vector<double>& D,
    const vector<double>& E,
    const vector<double>& b,
    vector<double>& results
)
{
    int N = D.size();
    vector<double> P(N);
    vector<double> Q(N);

    P[0] = -E[0]/D[0];
    Q[0] = b[0]/D[0];
    for(int i = 1; i < N; i++){
        P[i] = -E[i] / (D[i]+C[i]*P[i-1]);
        Q[i] = (b[i] - C[i]*Q[i-1])/(D[i]+C[i]*P[i-1]);
    }

    results[N-1] = Q[N-1];
    for (int i = N-2; i >= 0; i--){
        results[i] = P[i]*results[i+1] + Q[i];
    }
}