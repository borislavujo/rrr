function vsl = logSumDiff(Trms)
n = size(Trms,1);
Trms = sortrows(Trms);
vsl = cistLogSumDiff(Trms);
if (vsl(1)-Trms(n,1)<1e6)
  Trms = flipud(Trms);
  vsl = cistLogSumDiff(Trms);
endif
endfunction


function vsl = cistLogSumDiff(Trms)
n = size(Trms,1);
vsl = Trms(1,:);
for i=2:n
  if (vsl(2)==Trms(i,2))
    vsl(1) = logSumExp([vsl(1);Trms(i,1)]);
  elseif (vsl(1)>Trms(i,1))
    vsl(1) = logDiffExp(vsl(1),Trms(i,1));
  elseif (vsl(1)<Trms(i,1))
    vsl(1) = logDiffExp(Trms(i,1),vsl(1));
    vsl(2) = -vsl(2);
  else
    vsl = [-1e9, 1];
  endif
endfor
endfunction