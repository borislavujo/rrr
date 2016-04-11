function vln = multLogMat3(vl,vlp)
% returns logs of elements T(1,2), T(1,3), T(2,3) of a transition matrix
% that is the square of the transition matrix given on input (vl)
l12 = vl(1); l13=vl(2); l23=vl(3);
lp1 = vlp(1); lp2 = vlp(2); lp3 = vlp(3);
vln = vl;
l21 = l12+lp2-lp1; % dependent non-digonal terms of the transition matrix
l31 = l13+lp3-lp1;
l32 = l23+lp3-lp2;
l11 = real(subFrom1(logSumExp([l21;l31]))); % diagonal terms
l22 = real(subFrom1(logSumExp([l12;l32])));
l33 = real(subFrom1(logSumExp([l13;l23]))); % why would any of these become imaginary?
vln(1) = logSumExp([l11+l12;l12+l22;l13+l32]);
vln(2) = logSumExp([l11+l13;l12+l23;l13+l33]);
vln(3) = logSumExp([l21+l13;l22+l23;l23+l33]);
endfunction
