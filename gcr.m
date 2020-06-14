function [x, converged, iter_cnt, res_norms] = gmres(A, b, res_tol, m)
  n = size(A, 1);
  x = zeros(n, 1);  
  r = b - A * x;
  p=zeros(n,m);
  beta=zeros(m,m)
	residual = norm(r);
	stop_res = residual * res_tol;
  out_iter = 0;
	iter_cnt = 1;
	res_norms = zeros(max_iter * restart + 1, 1);
	res_norms(iter_cnt) = residual;
  converged = 0;
  p(:,1)=r
  pom=zeros(n,1);
  for j=1:m
    alpha=(r'*A*p(:,j))./((A*p(:,j))'*(A*p(:,j)));
    x=x+alpha*p(:,j);
    r=r-alpha*A*p(:,j);
    for i=1:j
      beta=-((A*r)'*A*p(:,i))./((A*p(:,i))'*A*p(:,i))*p(:,i);
      pom=pom+beta
    endfor
    p(:,j+1)=r+pom
    pom=zeros(n,1);
  endfor
