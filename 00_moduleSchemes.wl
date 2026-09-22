#!/usr/bin/env wolframscript

GetFiniteDifference[xPoints_List, fPoints_List, derivOrder_Integer] := Module[{n, c, p, eqns, sol, scheme},
  
  (*getting the polynomial order*)
  n = Length[xPoints]; 
  (*number of coefficients needed for a (n-1)^th degree polynomial*)
  
  (*generating a list of coefficients C0, C1, C2,..Cn-1 *)
  coeffs = Array[c, n, 0];
  
  (*defining the n-1^th degree polynomial: C0 + C1*x + C2*x^2 + .. + Cn-1*x^n-1*)
  p[x_] = Sum[coeffs[[i + 1]]*x^i, {i, 0, n - 1}];
  
  (*setting up the system of equations*)
  eqns = Table[p[xPoints[[i]]] == fPoints[[i]], {i, 1, n}];
  
  (*solving for the coefficients*)
  sol = Solve[eqns, coeffs];
  
  (*extracting the required coefficients with numerical prefactor 
  based on the required derivative*)
  scheme = (derivOrder!)*coeffs[[derivOrder + 1]] /. sol[[1]];
  
  (*simplifying the result*)
  Simplify[scheme]]