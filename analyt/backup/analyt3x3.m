function [lkab,lkba] = analyt3x3(Kl,vpl)
% calculates analytically rates between states (1,2) and (3)
vpl = vpl - logSumExp(vpl);
% 1. make state a the one of {1,2} having smaller rxn time with (3)
if (Kl(1,3)<Kl(2,3))
  P = [0,1,0;1,0,0;0,0,1];
  vpl = P*vpl
  Kl = P*Kl*P
endif
% 2. subtract from Kl inverse rxn time tac
ltac = logSumExp([Kl(1,3);Kl(3,1)]);
Kl = Kl - ltac;
% 3. assign pa,pb,tab,tbc
lpa = vpl(1);
lpb = vpl(2);
ltab = logSumExp([Kl(1,2);Kl(2,1)]);
ltbc = logSumExp([Kl(2,3);Kl(3,2)]);
% 4. calculate the terms for numerator
vpref = [zeros(6,1);log(2);log(3);0;log(3);zeros(5,1)];
vtab = [ones(4,1);zeros(7,1);ones(4,1)];
vtbc = [zeros(4,1);ones(11,1)];
vpa = zeros(15,1); vpa([3;7;10;12;14]) = 1; vpa([1;4;5;8]) = 2;    vpa([2;6]) = 3;
vpb = zeros(15,1); vpb([3;4;7;8;12]) = 1;   vpb([9;10;13;14]) = 2; vpb([11;15]) = 3;
vppm = ones(15,1); vppm([2;4;6;8;10;11;14;15]) = -1;
vhod = vpref + vtab*ltab + vtbc*ltbc + vpa*lpa + vpb*lpb;
Ujo = [vhod,vppm];
vnum = logSumDiff(Ujo)
% 5. calculate the terms for denominator
vpref = zeros(17,1); vpref([2;9;15]) = log(2);
vtab = [zeros(5,1);ones(7,1);zeros(5,1)];
vtbc = [zeros(12,1);ones(5,1)];
vpa = ones(17,1); vpa([1;4;8;11]) = 0; vpa([3;7;10;14;16]) = 2;
vpb = ones(17,1); vpb([6;7;13;14]) = 0; vpb([4;5;11;12;17]) = 2;
vppm = ones(17,1); vppm([2;4;7;9;11;14;15]) = -1;
vhod = vpref + vtab*ltab + vtbc*ltbc + vpa*lpa + vpb*lpb;
%format long
Ujo = [vhod,vppm]
vden = logSumDiff(Ujo)
% 6. calculate the rate (add what was subtracted)
lkab = ltac + vnum(1) - vden(1)
lkba = lkab + vpl(3) - logSumExp(vpl(1:2))
