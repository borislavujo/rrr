function ld = logDiffExp(a,b)
% returns log (exp a - exp b)
d = a-b;
ld = real(a + subFrom1(-d));
