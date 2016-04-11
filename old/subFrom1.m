function lp = subFrom1(l1,howFine)
% subtract a log number from 1
if (nargin<2) howFine = 10; endif
if (l1<-1e-5) 
  lp = real(log(1-exp(l1)));
else
  p = 0;
  for i=1:howFine
    p-=l1^i/factorial(i);
  endfor
  lp = real(log(p));
endif
endfunction
