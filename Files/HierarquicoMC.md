param ni > 0; 	# numero de clientes
param nj > 0; 	# numero de facilidades

param p > 0; 	# numero de facilidades disponiveis
param q > 0;	# numero de facilidades de segundo nivel disponiveis

set I := {0..ni-1};		# conjunto de clientes
set J := {0..nj-1};		# conjunto de facilidades

param f{I};			# demanda de cada cliente

param a{I,J} >= 0;	# matriz de atendimento da facilidade de nivel baixo
param b{I,J} >= 0;	# matriz de atendimento caso haja confronto de facilidades
param c{I,J} >= 0;	# matriz de atendimento da facilidade de nivel alto

var x{I}, binary;	# decisao do cliente ser ou nao atendido
var y{J}, binary;	# decisao da facilidade de nivel baixo ser ou nao utilizada
var z{J}, binary;	# decisao da facilidade de nivel alto ser ou nao utilizada


maximize fo : sum{i in I} f[i]*x[i];

s.t. r1{i in I}: sum{j in J} a[i,j] *y[j] + sum{j in J} b[i,j] * z[j] - x[i] >= 0;
s.t. r2{i in I}: sum{j in J} c[i,j]*z[j] - x[i] >= 0;
s.t. r3: sum{j in J} y[j] = p;
s.t. r4: sum{j in J} z[j] = q;