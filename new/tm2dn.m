function [D,N] = tm2dn(L,vpl)
n = size(L,1);
D = zeros(n);
N = ones(n);
for i=1:n
  for j=1:n
    if (L(i,j)>vpl(i))
      D(i,j) = subFrom1(vpl(i) - L(i,j)) + L(i,j) - vpl(i);
    else
      D(i,j) = subFrom1(L(i,j) - vpl(i));
      N(i,j) = -1;
    endif
  endfor
endfor
endfunction