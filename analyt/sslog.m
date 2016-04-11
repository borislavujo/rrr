function [kabl, kbal] = sslog(Kl, vpl)
% rate constnt from the steady state approximation
Kl(find(Kl==0))=-1e9;
vpl = vpl - logSumExp(vpl);

pl12 = logSumExp(vpl(1:2));
ppl1 = vpl(1)-pl12;
ppl2 = vpl(2)-pl12;
ppl12 = pl12;
ppl3 = vpl(3);

pl11 = 0;
pl12 = Kl(2,1)-logSumExp([Kl(3,2);Kl(1,2)]);
pl1 = logSumExp([pl11;pl12]);
fl11 = Kl(3,1);
fl12 = Kl(3,2)+pl12;
fl1 = logSumExp([fl11;fl12]);
tel1 = pl1 - fl1; % log exit time from state 1

pl21 = Kl(1,2) - logSumExp([Kl(3,1);Kl(2,1)]);
pl22 = 0;
pl2 = logSumExp([pl21;pl22]);
fl21 = Kl(3,1) + pl21;
fl22 = Kl(3,2);
fl2 = logSumExp([fl21;fl22]);
tel2 = pl2 - fl2; % log exit time from state 2

kbal = -max(tel1,tel2); % 1/longer time is the lowest possible rate constant from (1+2) to 3
kabl = kbal+ppl12-ppl3; % detailed balance is satisfied
%kabtst = tstlog(Kl,vpl);
%tstss = exp(kabtst-kabl)