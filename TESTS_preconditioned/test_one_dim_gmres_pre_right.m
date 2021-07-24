addpath('C:\Users\rysza\Desktop\PMPMPMPM\implementation\GMRES\preconditioning')
addpath('C:\Users\rysza\Desktop\PMPMPMPM\implementation\GCR\preconditioning')
addpath('C:\Users\rysza\Desktop\PMPMPMPM\implementation\ORTHODIR\preconditioning')
addpath('C:\Users\rysza\Desktop\PMPMPMPM\implementation\ORTHOMIN\preconditioning')
addpath('C:\Users\rysza\Desktop\PMPMPMPM\implementation\Matrices')
addpath('C:\Users\rysza\Desktop\PMPMPMPM\implementation\Conditioners')

b=1;
c=1;
d=0;
e=pi;
f=@(x)2.*sin(x)+cos(x);
N=151;
alpha=0;
beta=0;
stop_res=0.0001;
max_iter=300;
m_gmres=10;

omega=0.5;

[A,b]=one_dim(b,c,d,e,alpha,beta,f,N);
[M1]=jacobi(A);
[M2]=gauss_seidel(A);
[M3]=sor(A,omega);
[M4]=incompletelu(A);

[x1_right_M1, converged, ic_gmres_right_M1, rn_gmres_right_M1] = gmres_pre_right(A, b,m_gmres, max_iter, stop_res, M1);
[x1_right_M2, converged, ic_gmres_right_M2, rn_gmres_right_M2] = gmres_pre_right(A, b,m_gmres, max_iter, stop_res, M2);
[x1_right_M3, converged, ic_gmres_right_M3, rn_gmres_right_M3] = gmres_pre_right(A, b,m_gmres, max_iter, stop_res, M3);
[x1_right_M4, converged, ic_gmres_right_M4, rn_gmres_right_M4] = gmres_pre_right(A, b,m_gmres, max_iter, stop_res, M4);

semilogy(0 : ic_gmres_right_M1-1, rn_gmres_right_M1, 'r-'), hold on
semilogy(0 : ic_gmres_right_M2-1, rn_gmres_right_M2,  'm-'), hold on
semilogy(0 : ic_gmres_right_M3-1, rn_gmres_right_M3,  'g-'), hold on
semilogy(0 : ic_gmres_right_M4-1, rn_gmres_right_M4,  'c-'), hold on

xmin = 0;
xmax = max([ic_gmres_right_M1+1, ic_gmres_right_M2+1, ic_gmres_right_M3+1, ic_gmres_right_M4+1]) + 1;
ymin = min([min(rn_gmres_right_M1), min(rn_gmres_right_M2), min(rn_gmres_right_M3), min(rn_gmres_right_M4)]) * 0.8;
ymax = max([max(rn_gmres_right_M1), max(rn_gmres_right_M2), max(rn_gmres_right_M3), max(rn_gmres_right_M4)]) * 1.2;
axis([xmin xmax ymin ymax]);

xlabel('Liczba Iteracji'), ylabel('2 norma residuum'), hold on
legend('Jacobi', 'Gauss-Seidel', 'SOR(1/2)','ILU(0)'), hold on
title_str1 = 'GMRES z prawym preconditioningiem dla jednowymiarowego zadania';
title({title_str1}), hold off
saveas(gcf,'C:\Users\rysza\Desktop\PMPMPMPM\implementation\Plots\preconditioned\right_preconditioned_problem1.jpg')