function K = expmat(Kl)
n = size(Kl,1);
K = Kl;
for i=1:n
  K(i,find(Kl(i,:))) = exp(Kl(i,find(Kl(i,:))));
endfor
K = full(K - diag(sum(K)));