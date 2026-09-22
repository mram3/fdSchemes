#include "ExplicitSchemes.h"
#include <stdexcept>
#include<cmath>

double CD4A::firstDeriv(
    const Grid& mesh, int i
) const
{
    double h = mesh.x[i]-mesh.x[i-1];
    double phi1 = mesh.phi[i-2];
    double phi2 = mesh.phi[i-1];
    double phi4 = mesh.phi[i+1];
    double phi5 = mesh.phi[i+2];

    double a = phi1 + 8 * (phi4-phi2) - phi5;

    return a / (12 * h);
}

double CD4A::secondDeriv(
    const Grid& mesh, int i
) const
{
    double h = mesh.x[i]-mesh.x[i-1];
    double phi1 = mesh.phi[i-2];
    double phi2 = mesh.phi[i-1];
    double phi3 = mesh.phi[i];
    double phi4 = mesh.phi[i+1];
    double phi5 = mesh.phi[i+2];

    double a = - 30*phi3 + 16 * (phi2+phi4) - phi1 - phi5;

    return a / (12 * h * h);
}

double CD2A::firstDeriv(
    const Grid& mesh, int i
) const
{
    double h1 = mesh.x[i] - mesh.x[i-1];
    double h2 = mesh.x[i+1] - mesh.x[i];
    double phi1 = mesh.phi[i-1];
    double phi2 = mesh.phi[i];
    double phi3 = mesh.phi[i+1];

    double a = phi3*h1*h1 - phi1*h2*h2 + phi2*(h2*h2-h1*h1);
           a/= h1*h2*(h1+h2);
    return a;
}

double CD2A::secondDeriv(
    const Grid& mesh, int i
) const
{
    double h1 = mesh.x[i]   - mesh.x[i-1]; 
    double h2 = mesh.x[i+1] - mesh.x[i];
    double h3 = mesh.x[i+2] - mesh.x[i+1];
    double phi1 = mesh.phi[i-1];
    double phi2 = mesh.phi[i];
    double phi3 = mesh.phi[i+1];
    double phi4 = mesh.phi[i+2];

    double a = phi1*h2*h3*(2*h2*h2+3*h2*h3+h3*h3)
              +phi2*(h1+h2)*h3*(h1*h1-h1*h2-2*h2*h2-3*h2*h3-h3*h3)
              +phi3*h1*(h2+h3)*(-h1*h1+pow(h2+h3, 2))
              +phi4*h1*h2*(h1*h1-h2*h2);
    
           a /= h1*h2*(h1+h2)*h3*(h2+h3)*(h1+h2+h3);

    return 2 * a;
}

double CD1A::firstDeriv(
    const Grid& mesh, int i
) const
{
    throw std::runtime_error("Error: First order accurate CDS does not exist");
    return 0.0;
}
double CD1A::secondDeriv(
    const Grid& mesh, int i
) const
{
    double h1 = mesh.x[i] - mesh.x[i-1];
    double h2 = mesh.x[i+1] - mesh.x[i];
    double phi1 = mesh.phi[i-1];
    double phi2 = mesh.phi[i];
    double phi3 = mesh.phi[i+1];

    double a =  phi3*h1 + phi1*h2 - phi2*(h1+h2);
           
           a/= h1*h2*(h1+h2);
    
    return 2 * a;
}

double FD4A::firstDeriv(
    const Grid& mesh, int i
) const
{
    double h = mesh.x[i+1]-mesh.x[i];
    double phi1 = mesh.phi[i];
    double phi2 = mesh.phi[i+1];
    double phi3 = mesh.phi[i+2];
    double phi4 = mesh.phi[i+3];
    double phi5 = mesh.phi[i+4];

    double a = - 25*phi1 + 48*phi2 - 36*phi3 + 16*phi4 - 3*phi5;

    return a / (12 * h);
}

double FD4A::secondDeriv(
    const Grid& mesh, int i
) const
{
    double h = mesh.x[i+1] - mesh.x[i];
    double phi1 = mesh.phi[i];
    double phi2 = mesh.phi[i+1];
    double phi3 = mesh.phi[i+2];
    double phi4 = mesh.phi[i+3];
    double phi5 = mesh.phi[i+4];
    double phi6 = mesh.phi[i+5];

    double a = 45*phi1 - 154*phi2 + 214*phi3 - 156*phi4 + 61*phi5 - 10*phi6;

    return a / (12 * h * h);
}

double FD2A::firstDeriv(
    const Grid& mesh, int i
) const
{
    double h1 = mesh.x[i+1] - mesh.x[i];
    double h2 = mesh.x[i+2] - mesh.x[i+1];
    double phi1 = mesh.phi[i];
    double phi2 = mesh.phi[i+1];
    double phi3 = mesh.phi[i+2];

    double a = -phi3*h1*h1 + phi2*pow(h1+h2,2) - phi1*h2*(2*h1+h2);
           a/= h1*h2*(h1+h2);
    return a;
}

double FD2A::secondDeriv(
    const Grid& mesh, int i
)const
{
    double h1 = mesh.x[i+1] - mesh.x[i]; 
    double h2 = mesh.x[i+2] - mesh.x[i+1];
    double h3 = mesh.x[i+3] - mesh.x[i+2];
    double phi1 = mesh.phi[i];
    double phi2 = mesh.phi[i+1];
    double phi3 = mesh.phi[i+2];
    double phi4 = mesh.phi[i+3];

    double a = phi1*h2*h3*(h2+h3)*(3*h1+2*h2+h3)
              -phi2*(h1+h2)*h3*(2*h1*h1+4*h1*h2+2*h2*h2+3*h1*h3+3*h2*h3+h3*h3)
              +phi3*h1*(h2+h3)*(2*h1*h1+3*h1*(h2+h3)+pow(h2+h3, 2))
              -phi4*h1*h2*(2*h1*h1+h2*h2+3*h1*h2);
    
           a /= h1*h2*(h1+h2)*h3*(h2+h3)*(h1+h2+h3);

    return 2 * a;
}

double FD1A::firstDeriv(
    const Grid& mesh, int i
)const
{
    double h1 = mesh.x[i+1] - mesh.x[i];
    double phi1 = mesh.phi[i];
    double phi2 = mesh.phi[i+1];

    return (phi2-phi1)/h1;
}

double FD1A::secondDeriv(
    const Grid& mesh, int i
)const
{
    double h1 = mesh.x[i+1] - mesh.x[i];
    double h2 = mesh.x[i+2] - mesh.x[i+1];
    double phi1 = mesh.phi[i];
    double phi2 = mesh.phi[i+1];
    double phi3 = mesh.phi[i+2];

    double a =  phi3*h1 + phi1*h2 - phi2*(h1+h2);
           
           a/= h1*h2*(h1+h2);
    
    return 2 * a;
}

double BD4A::firstDeriv(
    const Grid& mesh, int i
) const
{
    double h = mesh.x[i]-mesh.x[i-1];
    double phi1 = mesh.phi[i];
    double phi2 = mesh.phi[i-1];
    double phi3 = mesh.phi[i-2];
    double phi4 = mesh.phi[i-3];
    double phi5 = mesh.phi[i-4];

    double a = 25*phi1 - 48*phi2 + 36*phi3 - 16*phi4 + 3*phi5;

    return a / (12 * h);
}

double BD4A::secondDeriv(
    const Grid& mesh, int i
) const
{
    double h = mesh.x[i] - mesh.x[i-1];
    double phi1 = mesh.phi[i];
    double phi2 = mesh.phi[i-1];
    double phi3 = mesh.phi[i-2];
    double phi4 = mesh.phi[i-3];
    double phi5 = mesh.phi[i-4];
    double phi6 = mesh.phi[i-5];

    double a = 45*phi1 - 154*phi2 + 214*phi3 - 156*phi4 + 61*phi5 - 10*phi6;

    return a / (12 * h * h);
}

double BD2A::firstDeriv(
    const Grid& mesh, int i
) const
{
    double h1 = mesh.x[i] - mesh.x[i-1];
    double h2 = mesh.x[i-1] - mesh.x[i-2];
    double phi1 = mesh.phi[i];
    double phi2 = mesh.phi[i-1];
    double phi3 = mesh.phi[i-2];

    double a = phi3*h1*h1 - phi2*pow(h1+h2,2) + phi1*h2*(2*h1+h2);
           a/= h1*h2*(h1+h2);
    return a;
}

double BD2A::secondDeriv(
   const Grid& mesh, int i
)const
{
    double h1 = mesh.x[i]   - mesh.x[i-1]; 
    double h2 = mesh.x[i-1] - mesh.x[i-2];
    double h3 = mesh.x[i-2] - mesh.x[i-3];
    double phi1 = mesh.phi[i];
    double phi2 = mesh.phi[i-1];
    double phi3 = mesh.phi[i-2];
    double phi4 = mesh.phi[i-3];

    double a = phi1*h2*h3*(h2+h3)*(3*h1+2*h2+h3)
              -phi2*(h1+h2)*h3*(2*h1*h1+4*h1*h2+2*h2*h2+3*h1*h3+3*h2*h3+h3*h3)
              +phi3*h1*(h2+h3)*(2*h1*h1+3*h1*(h2+h3)+pow(h2+h3, 2))
              -phi4*h1*h2*(2*h1*h1+h2*h2+3*h1*h2);
    
           a /= h1*h2*(h1+h2)*h3*(h2+h3)*(h1+h2+h3);

    return 2 * a;
}

double BD1A::firstDeriv(
    const Grid& mesh, int i
)const
{
    double h1 = mesh.x[i] - mesh.x[i-1];
    double phi1 = mesh.phi[i];
    double phi2 = mesh.phi[i-1];

    return (phi1 - phi2)/h1;
}

double BD1A::secondDeriv(
    const Grid& mesh, int i
)const
{
    double h1 = mesh.x[i] - mesh.x[i-1];
    double h2 = mesh.x[i-1] - mesh.x[i-2];
    double phi1 = mesh.phi[i];
    double phi2 = mesh.phi[i-1];
    double phi3 = mesh.phi[i-2];

    double a =  phi3*h1 + phi1*h2 - phi2*(h1+h2);
           
           a/= h1*h2*(h1+h2);
    
    return 2 * a;
}