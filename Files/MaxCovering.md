param ni > 0; 	# n�mero de clientes
param nj > 0; 	# n�mero de facilidades
param p > 0; 	# n�mero de facilidades dispon�veis

set I := {0..ni-1};		# conjunto de clientes
set J := {0..nj-1};		# conjunto de facilidades

param w{I};			# peso de cada cliente
set a{I} within J;	# matriz de facilidades que conseguem atender cada cliente

var x{J}, binary; 	# vari�vel que decide instalar a facilidade
var y{I}, binary;	# vari�vel que decide cobrir um cliente

maximize fo : sum{i in I} w[i]*y[i];		# maximizar os clientes cobertos

s.t. r1{i in I}: sum{j in a[i]} x[j] >= y[i];
s.t. r2: sum{j in J} x[j] = p;

end;