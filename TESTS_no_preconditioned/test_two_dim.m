addpath('C:\Users\rysza\Desktop\PMPMPMPM\implementation\GMRES')
addpath('C:\Users\rysza\Desktop\PMPMPMPM\implementation\GCR')
addpath('C:\Users\rysza\Desktop\PMPMPMPM\implementation\ORTHODIR')
addpath('C:\Users\rysza\Desktop\PMPMPMPM\implementation\ORTHOMIN')
addpath('C:\Users\rysza\Desktop\PMPMPMPM\implementation\Matrices')

b1 = 0;
b2 = 0;
c = 0;
d = 0;
e = 1;
f = 0;
g = 1;
N=11;
fun = @(x,y)0;
h = @(x,y)sin(x).*cos(x);

[A,b]=two_dim(b1,b2,c,d,e,f,g,N,fun,h);


stop_res=0.0001;
max_iter=300;
m_gmres=10;
m=5;
k=m-3;

[x, converged, ic_gmres, rn_gmres] = gmres(A, b,m_gmres, max_iter, stop_res);
[x, converged,ic_gcr,rn_gcr] = gcr(A, b, m, max_iter,stop_res);
[x, converged,ic_orthomin,rn_orthomin] = orthomin(A, b, m, max_iter,stop_res,k);
[x, converged,ic_orthodir,rn_orthodir] = orthodir(A, b, m, max_iter,stop_res);