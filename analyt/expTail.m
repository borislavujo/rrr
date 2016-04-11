function ltau = expTail(vTl,vKl,vlp,ltau)
% integrates an estimate for the tail of the relaxation
l33 = subFrom1(logSumExp([vTl(2),vTl(3)]));
f13 = l33 + vKl(1);
f23 = l33 + vKl(2);
f31 = vTl(2) + vKl(3);
f32 = vTl(3) + vKl(4);
lflux = logDiffExp(logSumExp([f13,f23]),logSumExp([f31,f32]));
ldtau = 2 * logDiffExp(l33,vlp(3)) - lflux;
ltau = logSumExp(ltau,ldtau);
