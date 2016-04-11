function [L3,D3,N3] = multLogMat(vpl,L1,D1,N1,L2,D2,N2)
% returns transition matrix represented in the best formalism T_ij = p_i^eq (1 +/- e^D_ij)
% L_ij = log (T_ij)
% h decides whether L or D is calculated precisely
h = -3;
if (nargin<5) L2=L1; D2=D1; N2=N1; endif
n = size(vpl,1);
L3 = zeros(n); D3 = zeros(n); N3 = ones(n);
% first do L
for i=1:n
  for j=1:n
    L3(i,j) = logSumExp(L1(i,:)'+L2(:,j));
  endfor
endfor
for j=1:n
  L3(:,j) = L3(:,j) - logSumExp(L3(:,j)); % normalise L
endfor
% now D
for i=1:n
  for j=1:n
    if (L3(i,j)-vpl>h)
      Trms = zeros(n*3,2); % the second column is 1 if the term is to be subtracted
      for k=1:n
	inow = (k-1)*3+1;
	Trms(inow:inow+2,1) = vpl(k)*[1;1;1] + D1(i,k)*[1;0;1] + D2(k,j)*[0;1;1];
	Trms(inow:inow+2,2) = [N1(i,k);N2(k,j);N1(i,k)*N2(k,j)];
      endfor
      vddl = logSumDiff(Trms);
      D3(i,j) = vddl(1);
      N3(i,j) = vddl(2);
      if (N3(i,j)==1) L3(i,j) = vpl(i) + logSumExp([0;D3(i,j)]); else L3(i,j) = vpl(i) + subFrom1(D3(i,j)); endif
    else
      if (L3(i,j)>vpl(i)) N3(i,j) = 1; else N3(i,j) = -1; endif
      D3(i,j) = log(N3(i,j)*(exp(L3(i,j)-vpl(i))-1));
    endif
  endfor
endfor
% normalise
D3 = (D3+D3')/2;
for j=1:n
  Trms = [D3(:,j)+vpl, N3(:,j)];
  vcol = logSumDiff(Trms);
  for i=1:n
    vpozmene = logSumDiff([min(vcol(1),D3(i,j)-log(2)),-vcol(2);D3(i,j),N3(i,j)]);
    D3(i,j) = vpozmene(1);
  endfor
  ujo = logSumExp(L3(:,n));
  L3(:,n) = L3(:,n) - ujo;
endfor
endfunction

