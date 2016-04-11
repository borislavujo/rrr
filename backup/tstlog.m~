function [kABl, kBAl] = tstlog(Kl,vpl)
% calculates TST rates between states (1+2) and 3 
% 3x3 log rate matrix on input
kABl = logSumExp(Kl(1:2,3));
kBAl = kABl + vpl(3) - logSumExp(vpl(1:2));
endfunction
