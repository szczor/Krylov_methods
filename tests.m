A = [ 2 1  1;
     -1 1 -1;
      1 2  3];
b = [2; 3; -10];

res_tol  = 1e-9;
max_iter=10;
restart=2;
useHH=1

[x, converged, iter_cnt, res_norms] = gmres(A, b, res_tol, max_iter,restart,useHH)