addpath('C:\Users\rysza\Desktop\PMPMPMPM\implementation\GMRES\no_preconditioning')
addpath('C:\Users\rysza\Desktop\PMPMPMPM\implementation\GCR\no_preconditioning')
addpath('C:\Users\rysza\Desktop\PMPMPMPM\implementation\ORTHODIR\no_preconditioning')
addpath('C:\Users\rysza\Desktop\PMPMPMPM\implementation\ORTHOMIN\no_preconditioning')
addpath('C:\Users\rysza\Desktop\PMPMPMPM\implementation\Matrices')

b=1;
c=1;
d=0;
e=pi;
f=@(x)2.*sin(x)+cos(x);
N=101;
alpha=0;
beta=0;
stop_res=0.0001;
max_iter=300;
m_gmres=10;
m=5;
k=m-3;

[A,be]=one_dim(b,c,d,e,alpha,beta,f,N);
[x, converged, ic_gmres, rn_gmres] = gmres(A, be,m_gmres, max_iter, stop_res);
[x, converged,ic_gcr,rn_gcr] = gcr(A, be, m, max_iter,stop_res);
[x, converged,ic_orthomin,rn_orthomin] = orthomin(A, be, m, max_iter,stop_res,k);
[x, converged,ic_orthodir,rn_orthodir] = orthodir(A, be, m, max_iter,stop_res);

semilogy(0 : ic_gmres-1, rn_gmres, 'r-'), hold on
semilogy(0 : ic_gcr-1,  rn_gcr,  'y-'), hold on
semilogy(0 : ic_orthomin-1,  rn_orthomin,  'g-'), hold on
semilogy(0 : ic_orthodir-1,  rn_orthodir,  'c-'), hold on


xmin = 0;
xmax = max([ic_gmres+1, ic_gcr+1, ic_orthomin+1, ic_orthodir+1]) + 1;
ymin = min([min(rn_gmres), min(rn_gcr), min(rn_orthomin), min(rn_orthodir)]) * 0.8;
ymax = max([max(rn_gmres), max(rn_gcr), max(rn_orthomin), max(rn_orthodir)]) * 1.2;
axis([xmin xmax ymin ymax]);


xlabel('Liczba Iteracji'), ylabel('2 norma residuum'), hold on
legend('GMRES', 'GCR', 'ORTHOMIN', 'ORTHODIR'), hold on
title_str1 = 'Algorytmy bez preconditiningu dla jednowymiarowego zadania';
title({title_str1}), hold off

saveas(gcf,'C:\Users\rysza\Desktop\PMPMPMPM\implementation\Plots\not_preconditioned\not_preconditioned_problem1_105.jpg')