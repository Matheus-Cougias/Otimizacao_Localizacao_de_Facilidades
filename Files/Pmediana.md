param ni > 0;
param nj > 0;
param p > 0;

set I := {0..ni-1};
set J := {0..nj-1};

param c{I,J} >= 0;

var y{J}, binary;
var x{I,J} >= 0 <=1;

minimize fo : sum {i in I, j in J} c[i,j] * x[i,j];

s.t. r1{i in I}: sum{j in J} x[i,j] = 1;
s.t. r2{i in I, j in J}: x[i,j] <= y[j];
s.t. r3: sum{j in J} y[j] = p;