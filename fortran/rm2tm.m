function Tl = rm2tm(Kl, ldt)
% produces a log transition matrix Tlm for log lag time ldt with the table Tm1l of log(1-term) from a rate matrix Kl
% zeros are not replaced by -1e9
n = size(Kl,1);
Kl([1:n+1:n^2]) = -1e9;% change diagonal rates to 9                                                                              
if (ldt>-max(max(Kl))) error('too large step size'); endif
Tl = Kl + ldt;
for i=1:n
  mabyt0 = logSumExp(Tl(:,i));
  Tl(i,i) = subFrom1(mabyt0);% normalise, so that sum_i(T_ji) = 1
endfor
[vmaxes,vmi] = max(Tl,[],2);
endfunction