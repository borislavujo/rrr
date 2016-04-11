function [kAB, kBA] = rxnK(K,va,vpe,howFine,startdt,ratioEnd)
% calculates the rate constants between observable states A and B defined by va
% using rate matrix K between microstates with equilibrium populations vpe
% the cost is roughly twice the cost of only exponential integration
if (nargin<4) howFine = 4; endif % how fine is the integration
if (nargin<5) startdt = 1.5e-2; endif % the initial time step is startdt x (time of the fastest process in the network)
if (nargin<6) ratioEnd = 1e-3; endif % criterion for stopping the propagation: 1-ratioEnd relaxed for all groupings
n = size(K,1);
nFine = 2^howFine;
vpe = vpe/sum(vpe); % normalise equilibrium populations
K = symmetriseMat(K,vpe); % symmetrise the rate matrix
K = K - diag(sum(K)); % just to ensure the rate matrix is correct
p1 = sum(vpe.*va); p2 = 1-p1; % equilibrium population of states 1 and 2
vp0 = vpe.*va/p1; % p is 1 at the beginning
dt0 = startdt/max(abs(diag(K))); % quite arbitrary compromise between efficiency and accuracy
t = dt0/nFine;
G = (eye(n)+K*t); F = G; % transition matrix at dt0
for iMulti = 1:howFine
  F = symmetriseMat(F,vpe);
  F = F*F;
endfor
p = ((va'*F*vp0)-p1)/p2
tau = dt0*(1+p)/2; t = dt0;
% relaxation
while (p>ratioEnd)
  [dtau,p] = integRxn(t,nFine,F,G,va,vp0,p1);
  tau += dtau;
  t = t*2; % t -> 2*t
  F = F*F; % t -> 2*t
  G = G*G; % t -> 2*t
  G = symmetriseMat(G,vpe);
  G = G + diag(ones(1,n)-sum(G,1)); % sum of elements in each column must be exactly 1
  F = symmetriseMat(F,vpe);
  F = F + diag(ones(1,n)-sum(F,1)); % sum of elements in each column must be exactly 1
endwhile

vdp = K*(F*vp0); dp = sum(vdp.*va)./p2; k = -dp./p
tau = (tau + p./k)
kAB = p2/tau; kBA = p1/tau; % add the tail
endfunction

function T = symmetriseMat(T,vpe)
  n = size(vpe,1);
  G = diag(ones(n,1)./vpe)*T;
  T = diag(vpe)*(G+G')/2;
  T(abs(T) < sqrt(realmin)) = 0; % avoid complicated algebraic manipulations with numbers of order 1e-154 and less
endfunction

function [tau,p] = integRxn(t,nFine,F,G,va,vpe,p1)
vt = linspace(t,2*t,nFine+1);
vp = zeros(1,nFine+1);
vpe = F*vpe;
vp(1) = (va'*vpe-p1)/(1-p1);
for i = 1:nFine
  vpe = G*vpe;
  vp(i+1) = (va'*vpe-p1)/(1-p1);
endfor
tau = trapz(vt,vp);
p = vp(nFine+1);
endfunction