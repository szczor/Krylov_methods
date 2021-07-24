function [M]=sor(A,omega)
  x = diag(A);
  L = tril(A,-1);
  M = diag(x);
  M=M+omega.*L;
endfunction