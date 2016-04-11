function [irow,icol,kmax] = maxNZ(K);
% reutrns indices of the largest non-zero term in a square matrix
n = size(K,1);
irow = 0; icol = 0; kmax = -1e99;
for i=1:n
  vi = find(K(:,i));
  if (size(vi,1))
    vrow = K(vi,i);
    [kmaxn,irowi] = max(vrow);
    if (kmaxn>kmax)
      icol = i;
      irow = vi(irowi);
      kmax = full(kmaxn);
    endif
  endif
endfor
