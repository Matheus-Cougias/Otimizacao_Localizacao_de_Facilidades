param ni > 0; 	# n�mero de clientes
param nj > 0; 	# n�mero de facilidades
param p > 0; 	# n�mero de facilidades dispon�veis

set I := {0..ni-1};		# conjunto de clientes
set J := {0..nj-1};		# conjunto de facilidades

param w{I};			# peso de cada cliente
set a{I} within J;	# matriz de facilidades que conseguem atender cada cliente

var x{J}, binary; 	# vari�vel que decide instalar a facilidade
var y{I}, binary;	# vari�vel que decide cobrir um cliente

minimize fo : sum{i in I} w[i]*y[i];		# minimiza os clientes n�o cobertos
s.t. r1{i in I}: y[i] + sum{j in a[i]} x[j] >= 1;
s.t. r2: sum{j in J} x[j] = p;		# a quantidade total de facilidades instaladas devem ser menor que o limitador P

end;