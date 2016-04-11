function [Kl, vpl] = rlm(n,m1,m2)
% returns nxn log rate matrix, diagonal is 0
if (nargin<3) m2 = 20; endif % magnitude of the equil pop distribution
if (nargin<2) m1 = 20; endif % magnitude of the rxn time distribution
vpl = m2*rand(n,1); vpl = vpl - logSumExp(vpl);
Taul = m1*rand(n);
Kl = zeros(n);
for i1=1:n-1
  for i2 = i1+1:n
    sp = logSumExp([vpl(i1);vpl(i2)]);
    Kl(i1,i2) = vpl(i1) - sp - Taul(i1,i2);
    Kl(i2,i1) = vpl(i2) - sp - Taul(i1,i2);
  endfor
endfor
endfunction