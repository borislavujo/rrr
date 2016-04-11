function lsum = logSumExp(vl)
% returns logarithm of the sum of exponentials of elements of vl
[m,n] = size(vl);
if and(n<2,m<2) 
  lsum = vl; 
else 
  if (m==1) vl = vl'; m = size(vl,1); endif
  vl = sort(vl); % ascending
  lsum = vl(1);
  for i=2:m
    lnow = vl(i);
    lsum = lnow + log(1+exp(lsum-lnow));
  endfor
endif
endfunction