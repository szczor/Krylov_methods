function [V,H,beta] = arnoldi_pre_right(A,m,r_0, M)
% Zmodyfikowany G-S
  beta = norm(r_0, 2);
  N = size(A, 1);
  V = zeros(N,m+1); 
  H = zeros(m + 1,m);
  V(:, 1) = r_0 / beta;
  for j = 1 : m
      u = M \ V(:,j);
      w = A*u;
      for i = 1 : j
          H(i,j) = w' * V(:,i);
          w = w - H(i,j) * V(:,i);
      end
      H(j+1,j) = norm(w,2);
      if (H(j+1,j) == 0)
        break;
      end
      V(:,j+1) = w / H(j+1,j);
  end
       
endfunction