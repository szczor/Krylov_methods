function[A,b]=random(N)
  A = sprandn(N,N,0.5) + 12*speye(N);
  b = rand(N,1);