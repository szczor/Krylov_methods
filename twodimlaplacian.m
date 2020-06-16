#{
Ten program tworzy macierz A i be dla dwuwymiarowego zadania równania ró¿niczkowego
-laplacianu+b^T\grad u +cu=f
gdzie x_1 jest w przedziale [d,e] oraz x_2 w przedziale [f,g] 
#}
function[A,be]=twodimlaplacian(b,c,d,e,f,g,N,fun,h)

h1=(e-d)./N;
h2=(g-f)./N;
u=zeros(N+1,N+1);
for l=1:N
  u(1,l)=h(d,f+l*h2);
  u(l,N+1)=h(d+l*h1,g);
  u(N+1,l)=h(e,f+l*h2);
  u(l,1)=h(d+l*h1,f);
end
be=zeros(N-1,N-1);
for i=1:(N-1)
  for j=1:(N-1)
    x1=d+i*h1;
    x2=f+j*h2;
    be(i,j)=fun(x1,x2);
    if(i==1)
      be(i,j)=be(i,j)+h(d,f+j*h2)./(h2.^2)
    elseif(j==1)
      be(i,j)=be(i,j)+h(d+i.*h1,g)./(h1.^2)
    elseif(i==N-1)
      be(i,j)=be(i,j)+h(e,f+j*h2)./(h2.^2)
    elseif(j==N-1)
      be(i,j)=be(i,j)+h(d+i.*h1,f)./(h1.^2)
    end
  end
end
A=zeros(N-1,N-1);
end