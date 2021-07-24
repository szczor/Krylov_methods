function [M]=incompletelu(A)
  N=size(A,1);
  for k=1 : (N-1)
    for i=(k+1) : N
      if  A(i,k) ~= 0
        A(i,k)=A(i,k)./A(k,k);
      end
    end
    for j=(k+1) : N
      for i=(k+1) : N
        if A(i,j)~= 0
          A(i,j)=A(i,j)-A(i,k).*A(k,j);
        end
      end
    end
  end
  L=eye(N);
  L= L+tril(A,-1);
  U= triu(A);
  M=L*U;
endfunction