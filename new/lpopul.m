function x = lpopul(L,D,N,vpl,va)
% calculates numerator of progress of the relaxation (current population in state A - the equilibrium population)
% first try simply Tl
vi = find(va);
n = size(vi,1);
vpl = vpl(vi); L = L(vi,vi); D = D(vi,vi); N = N(vi,vi);
pl1 = logSumExp(vpl); pl2 = subFrom1(pl1);
vplnow = zeros(n,1);
for i=1:n
  vplnow(i) = logSumExp(L(i,:)'+vpl-pl1);
endfor
pl = logSumExp(vplnow);
x = logDiffExp(pl,pl1)-pl2
if (x>-3)
  return;
endif

From1 = zeros(n^2,2);
for i=1:n
  for j=1:n
    k = (i-1)*n+j;
    From1(k,:) = [vpl(i) + vpl(j) + D(i,j), N(i,j)];
  endfor
endfor
vojoj = logSumDiff(From1);
x = vojoj(1)-pl1-pl2
endfunction
