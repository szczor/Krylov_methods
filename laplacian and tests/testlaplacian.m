
#{
b=1
c=1
f=@(x)sin(x)
d=0
e=pi
alpha=0
beta=0
N=10
[A,b]=laplacian(b,c,d,e,alpha,beta,f,N)
[x, converged, iter_cnt, res_norms] = gmres(A, b, res_tol, max_iter,restart,useHH)
[x, converged,iter_cnt,res_norms] = gcr(A, b, res_tol, m)
[x, converged,iter_cnt,res_norms] = orthomin(A, b, res_tol, m,k)
[x, converged,iter_cnt,res_norms] = orthodir(A, b, res_tol, m)
#}

h=@(x,y)sin(x).*cos(x);
fun=@(x,y)sin(x).*cos(x).*4;
c=0;
N=10;
d=0;
e=1;
f=0;
g=1;
b=0;
[A,be]=twodimlaplacian(b,c,d,e,f,g,N,fun,h)