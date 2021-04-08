param ni > 0; 	# numero de clientes
param nj > 0; 	# numero de facilidades

param p > 0; 	# numero de facilidades disponiveis
param q > 0;	# numero de facilidades de segundo nivel disponiveis

set I := {0..ni-1};		# conjunto de clientes
set J := {0..nj-1};		# conjunto de facilidades de baixo n�vel
set K := {0..nj-1};		# conjunto de facilidades de alto n�vel

param a{I,J} >= 0;      # simplificacao dos custos para facilidade de baixo nivel
param b{I,K} >= 0;      # simplificacao dos custos para facilidade de alto nivel atendendo no alto nivel
param c{I,K} >= 0;      # simplificacao dos custos para facilidade de alto nivel atendendo no baixo nivel


var y{J}, binary;	# decisao da facilidade de nivel baixo ser ou nao utilizada
var z{K}, binary;	# decisao da facilidade de nivel alto ser ou nao utilizada
var m{I,J} >= 0 <= 1;	# decisao de qual facilidade de baixo n�vel atender� qual cliente
var n{I,K} >= 0 <= 1;	# decis�o de qual facilidade de alto n�vel atender� qual cliente
var o{I,K} >= 0 <= 1;	# decisao de se a facilidade de alto n�vel instalada conseguir� atender no baixo n�vel


minimize fo: sum {i in I, j in J} a[i,j] * m[i,j] + sum{i in I, k in K} b[i,k] * n[i,k] + sum {i in I, k in K} c[i,k] * o[i,k];

s.t. r1 {i in I}: sum{j in J} m[i,j] + sum{k in K} o[i,k] = 1;
s.t. r2 {i in I}: sum{k in K} n[i,k] = 1;
s.t. r3 {i in I, j in J}: m[i,j] <= y[j];
s.t. r4 {i in I, k in K}: n[i,k] <= z[k];
s.t. r5 {i in I, k in K}: o[i,k] <= z[k];
s.t. r6: sum{j in J} y[j] = p;
s.t. r7: sum{k in K} z[k] = q;