function [x, converged,iter_cnt,res_norms] = gcr(A, b, res_tol, m)
  n = size(A, 1);
  x = zeros(n, 1);  
  r = b - A * x;
  p=zeros(n,m);
  beta=zeros(m,m);
	stop_res = norm(r) * res_tol;
  converged = 0;
  
  
  p(:,1)=r;
  pom=zeros(n,1);
  for j=1:m
    alpha=(r'*A*p(:,j))./((A*p(:,j))'*(A*p(:,j)));
    x=x+alpha*p(:,j);
    r=r-alpha*A*p(:,j);
    for i=1:j
      beta=-((A*r)'*A*p(:,i))./((A*p(:,i))'*A*p(:,i))*p(:,i);
      pom=pom+beta;
    endfor
    p(:,j+1)=r+pom;
    pom=zeros(n,1);
    res_norms(j) = norm(r,2);
		if (res_norms(j) < stop_res)
			break;
		end
  endfor
  if(res_norms(j)<stop_res) converged=1; end
  res_norms=res_norms(1:j)';
  iter_cnt=j;
end
