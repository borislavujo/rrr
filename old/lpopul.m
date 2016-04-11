function p = lpopul(vl,vlp)
% calculates progress of the relaxation
% state 3 relaxes to 1 and 2
lpn = logSumExp(vl(2:3)'); % popul of 1 and 2 now
lpe = logSumExp(vlp(1:2)); % equil popul of 1 and 2
dlp = lpn-lpe;
if (dlp<-30)
  p = 1;
else
  p = 1 - exp(dlp);
endif
endfunction
