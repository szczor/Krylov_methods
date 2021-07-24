#{
Rozwi¹¿ zadania najmniejszych kwadratów dla GMRES
#}
function [y] = leastsquare(H, beta)
	m = size(H, 2);
	b = zeros(m + 1, 1); b(1) = beta;
  y = H \ b;
endfunction
