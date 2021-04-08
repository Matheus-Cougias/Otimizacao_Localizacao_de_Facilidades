param ni > 0;  # number of clients
param nj > 0;  # number of facilities

set I := {0..ni-1};
set J := {0..nj-1};

param c{J};
set a{I} within J;

var x{J}, binary;

minimize of : sum{j in J} c[j] * x[j];

s.t. r1{i in I}: sum{j in a[i]} x[j] >= 1;

end; 
