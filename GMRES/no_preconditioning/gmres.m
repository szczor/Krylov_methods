function [x, converged,iter_cnt,res_norms]=gmres(A, b,m, max_iter,stop_res)
  n = size(A, 1);
  x = zeros(n, 1);   
  r=b-A*x;
  iter_cnt = 1;
  res_norms = zeros(max_iter,1);
  res_norms(1) = norm(r,2);
  converged = 0;
  while ((iter_cnt < max_iter) && (norm(r,2) > stop_res))
    [V,H,beta] = arnoldi(A,m,r);
    [y] = leastsquare(H,beta);
    z = V(:,1:m) * y(1:m);
    x = x + z;
    r = b - A*x;
    iter_cnt = iter_cnt+1;
    res_norms(iter_cnt) = norm(r,2);
  endwhile
  if(res_norms(iter_cnt) < stop_res) converged = 1; end
  res_norms = res_norms(1:iter_cnt)';  
endfunction
