function lsum = logSumExp(vl)
% returns logarithm of the sum of exponentials of elements of vl
[m,n] = size(vl);
if and(n<2,m<2) 
  lsum = vl; 
else 
  if (m==1) vl = vl'; m = size(vl,1); endif
  vl = sort(vl,"descend");
  ml = vl(1); % max term
  vl = vl - ml;
  s = 1;
  for i=2:m
    s += exp(vl(i));
  endfor
  lsum = real(log(s)+ml);
endif
endfunction