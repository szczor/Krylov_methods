#{
Ten program tworzy macierz A i be dla dwuwymiarowego zadania równania ró¿niczkowego
-laplacianu+b^T\grad u +cu=f
gdzie x_1 jest w przedziale [d,e] oraz x_2 w przedziale [f,g] 
#}
function[A,be]=two_dim(b_1,b_2,c,d,e,f,g,N,fun,h)

h1 = (e - d) ./ N;
h2 = (g - f) ./ N;
u = zeros(N+1,N+1);
for l = 1 : N
  u(1,l) = h(d,f+l*h2);
  u(l,N+1) = h(d+l*h1,g);
  u(N+1,l) = h(e,f+l*h2);
  u(l,1) = h(d+l*h1,f);
end
be = zeros(N-1,N-1);
for i = 1 : (N-1)
  x1 = d+i.*h1;
  for j = 1 : (N-1)
    x2 = f+j*h2;
    be(i,j) = fun(x1,x2);
    if i == 1
      be(i,j) = fun(x1,x2)+h(d,f+j.*h2)./(h1.^2);
    elseif i == N-1
      be(i,j) = fun(x1,x2)+h(e,f+j*h2)./(h1.^2)+b_1.*h(e,f+j*h2)./h1;
    end  
    if j == 1
      be(i,j) = fun(x1,x2)+h(d+i.*h1,g)./(h2.^2);
    elseif j == (N-1)
      be(i,j) = fun(x1,x2)+h(d+i.*h1,f)./(h2.^2)+b_2.*h(d+i.*h1,f)./h2;
    end
  end
end
be = be(:);
A = zeros((N-1)*(N-1),(N-1)*(N-1));
pom1 = (2./(h1.^2)+2./(h2.^2)-b_1./h1-b_2./h2);
pom2 = -1./(h1.^2)+b_1./h1;
pom3 = -1./(h2.^2)+b_2./h2;
pom4 = -1./(h2.^2);
pom5 = -1./(h1.^2);
for i = 1 : (N-1).*(N-1)
  for j=1 : (N-1).*(N-1)
    if(i == j)
      A(i,j) = pom1;
    elseif(j == i+1)
      A(i,j) = pom3;
    elseif(j == i-1)
      A(i,j) = pom4;
    elseif(j == i+N-1)
      A(i,j) = pom2;
    elseif(j == i-N+1)
      A(i,j) = pom5;
    else
      A(i,j) = 0;
    end
  end   
end