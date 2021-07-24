#{
this program creates a matrix A and vector b for a differential equation
 -u''(x)+bu'(x)+cu(x)=f(x), omega=[d,e] u(d)=alpha u(e)=beta
#}
function[A,be]=one_dim(b,c,d,e,alpha,beta,f,N)

h=(e-d)./N;
pom1=-1.*ones(1,N-1);
pom2=(b.*h-1).*ones(1,N-1);
pom3=(2-b.*h+c.*h.^2).*ones(1,N-1);
A=diag(pom3(1:N-1))+diag(pom2(2:N-1),1)+diag(pom1(2:N-1),-1);
A=A./(h.^2);
be=zeros(1,N-1);
x=zeros(1,N-1);
for i=1:(N-1)
  x(i)=d+i.*h;
  be(i)=f(x(i));
endfor

be(1)=be(1)+alpha./(h.^2);
be(N-1)=be(N-1)+((1- b .*h).*beta)./(h.^2);
be=be';
end