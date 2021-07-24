function [x, converged,iter_cnt,res_norms] = gcr(A, b, m, max_iter,stop_res)
  n = size(A, 1);
  x = zeros(n, 1);  
  r = b - A * x;
  p = zeros(n,m);
  beta = zeros(m,m);
  converged = 0;
  iter_cnt = 1;
  res_norms = zeros(max_iter,1);
  res_norms(1) = norm(r,2);
  
  while ((iter_cnt < max_iter) && (norm(r,2) > stop_res)) 
    p(:,1) = r;
    pom = zeros(n,1);
    for j = 1 : m
      alpha = (r' * A * p(:,j)) ./ ((A * p(:,j))' * (A * p(:,j)));
      x = x + alpha * p(:,j);
      r= r - alpha * A * p(:,j);
      for i = 1 : j
        beta = -((A * r)' * A * p(:,i))./((A * p(:,i))' * A * p(:,i)) * p(:,i);
        pom = pom + beta;
      endfor
      p(:,j+1) = r + pom;
    endfor
    r = b - A*x;
    iter_cnt = iter_cnt + 1;
    res_norms(iter_cnt) = norm(r,2);
  endwhile
  if(res_norms(iter_cnt) < stop_res) converged = 1; end
  res_norms = res_norms(1:iter_cnt)'; 
endfunction
