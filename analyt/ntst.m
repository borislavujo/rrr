function [kABl, kBAl] = ntst(Kl,vpl,via,vib)
if (nargin<3) via = [1;2]; vib = 3; endif
% calculates TST rates between states (1+2) and 3 
% 3x3 log rate matrix on input
na = size(via,1); nb = size(vib,1);
vklab = zeros(nb,1); vklba = zeros(na,1);
KltoA = Kl(via,vib);
for i=1:nb
  vipp = find(KltoA(:,i));
  if (size(vipp,1)>0)
    vklab(i) = logSumExp(KltoA(vipp,i));
  else
    vklab(i) = -1e9;
  endif
endfor
kABl = logSumExp(vklab+vpl(vib))-logSumExp(vpl(vib));
KltoB = Kl(vib,via);
for i=1:na
  vipp = find(KltoB(:,i));
  if (size(vipp,1)>0)
    vklba(i) = logSumExp(KltoB(vipp,i));
  else
    vklba(i) = -1e9;
  endif
endfor
kBAl = logSumExp(vklba+vpl(via))-logSumExp(vpl(via));
endfunction
