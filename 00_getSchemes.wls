#!/usr/bin/env wolframscript

SetDirectory[DirectoryName[$InputFileName]];
Get["00_moduleSchemes.wl"];

Print["==============================================="];
Print[" FINITE DIFFERENCE SCHEMES (using polynomials) "];
Print["===============================================\n"];


Print["-----------------------------------------------"];
Print["             FIRST DERIVATIVE                  "];
Print["-----------------------------------------------"];

Print["\n>>> 1st Order (Non-Uniform) <<<\n"];
Print["Forward  (FD): ", GetFiniteDifference[{0, h1}, {fi, fp1}, 1], "\n"];
Print["Backward (BD): ", GetFiniteDifference[{-hm1, 0}, {fm1, fi}, 1], "\n"];
(* There is no first order CD scheme *)

Print["\n>>> 2nd Order (Non-Uniform) <<<\n"];
Print["Central  (CD): ", GetFiniteDifference[{-hm1, 0, h1}, {fm1, fi, fp1}, 1], "\n"];
Print["Forward  (FD): ", GetFiniteDifference[{0, h1, h1+h2}, {fi, fp1, fp2}, 1], "\n"];
Print["Backward (BD): ", GetFiniteDifference[{-(hm1+hm2), -hm1, 0}, {fm2, fm1, fi}, 1], "\n"];

Print["\n>>> 4th Order (Uniform Stencil) <<<\n"];
Print["Central  (CD): ", GetFiniteDifference[{-2*h, -h, 0, h, 2*h}, {fm2, fm1, fi, fp1, fp2}, 1], "\n"];
Print["Forward  (FD): ", GetFiniteDifference[{0, h, 2*h, 3*h, 4*h}, {fi, fp1, fp2, fp3, fp4}, 1], "\n"];
Print["Backward (BD): ", GetFiniteDifference[{-4*h, -3*h, -2*h, -h, 0}, {fm4, fm3, fm2, fm1, fi}, 1], "\n"];


Print["\n\n-----------------------------------------------"];
Print["             SECOND DERIVATIVE                 "];
Print["-----------------------------------------------"];

Print["\n>>> 1st Order (Non-Uniform) <<<\n"];
(* 2nd derivative requires 3 points minimum just to exist (which yields 1st order on non-uniform grids) *)
Print["Central  (CD): ", GetFiniteDifference[{-hm1, 0, h1}, {fm1, fi, fp1}, 2], "\n"];
Print["Forward  (FD): ", GetFiniteDifference[{0, h1, h1+h2}, {fi, fp1, fp2}, 2], "\n"];
Print["Backward (BD): ", GetFiniteDifference[{-(hm1+hm2), -hm1, 0}, {fm2, fm1, fi}, 2], "\n"];

Print["\n>>> 2nd Order (Non-Uniform) <<<\n"];
(* To guarantee 2nd order accuracy on a non-uniform grid, we must add one extra point *)
Print["Central  (CD): ", GetFiniteDifference[{-hm1, 0, h1, h1+h2}, {fm1, fi, fp1, fp2}, 2], "\n"];
Print["Forward  (FD): ", GetFiniteDifference[{0, h1, h1+h2, h1+h2+h3}, {fi, fp1, fp2, fp3}, 2], "\n"];
Print["Backward (BD): ", GetFiniteDifference[{-(hm1+hm2+hm3), -(hm1+hm2), -hm1, 0}, {fm3, fm2, fm1, fi}, 2], "\n"];

Print["\n>>> 4th Order (Uniform Stencil) <<<\n"];
Print["Central  (CD): ", GetFiniteDifference[{-2*h, -h, 0, h, 2*h}, {fm2, fm1, fi, fp1, fp2}, 2], "\n"];
Print["Forward  (FD): ", GetFiniteDifference[{0, h, 2*h, 3*h, 4*h, 5*h}, {fi, fp1, fp2, fp3, fp4, fp5}, 2], "\n"];
Print["Backward (BD): ", GetFiniteDifference[{-5*h, -4*h, -3*h, -2*h, -h, 0}, {fm5, fm4, fm3, fm2, fm1, fi}, 2], "\n"];

Print["\n==============================================="];