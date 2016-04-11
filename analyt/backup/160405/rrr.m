function [vi,FEl,SSl,TSTl,NFl] = rrr(Kl,vpl,nDesiredStates,vObserv)
% grouping minima stops when either a rate constant threshold or a number of states is reached
% returns matrices without diagonal terms
% grouping is determined by largest of the lowest rate constants calculated using no-flux BC
% Kl ... log of the rate matrix
% vpl ... log of equil populs
% nDesiredStates ... grouping ends when this number of states is reached
% vObserv ... determines the states that should not be grouped together, as they belong to different observables
  n = size(Kl,1); nStates = n; % the original number of states
  if (nargin<4) vObserv = ones(n,1); endif % if anything can be grouped, make vObserv vector
  vObsInd = unique(vObserv); nObserv = size(vObsInd,1); % how many observable states there are
  if (nDesiredStates<nObserv) nDesiredStates = nObserv; endif % number of states at the end is at least the number of observables
  vi = [1:n]'; % index vector, always of size n
  Kl = Kl - diag(diag(Kl)); % delete all diagonal (negative) terms; Kl is now the decision matrix
  Kl = Kl -1e19*eye(n); % fill the diagonal, so it does not interfere with the searches for the largest rate constant
  SSl = Kl; % steady state matrix = the lower bound
  FEl = Kl; % fast equilibrium BC - the best upper bound
  NFl = Kl; % no-flux BC - should be between SSl and FEl
  TSTl = Kl; % TST matrix the worst upper bound
% avoiding grouping states from different observables
  if (nObserv>1)
    for i1Obs = 1:nObserv-1
      for i2Obs = i1Obs+1:nObserv
	vObs1 = find(vObserv==vObsInd(i1Obs))';
	vObs2 = find(vObserv==vObsInd(i2Obs))';
	for i1 = vObs1
	  for i2 = vObs2
	    if (Kl(i1,i2)~=0)
	      Kl(i1,i2) -= 1e9;
	      Kl(i2,i1) -= 1e9;
	    endif
	  endfor
	endfor
      endfor
    endfor
  endif
  [irow,icol,know] = maxNZ(min(Kl,Kl'));
  while (nStates>nDesiredStates)
    [Kl,vpln,vin]   = groupStates(Kl,vi,vpl,irow,icol,"nf"); % new decision matrix, assignment, and population vector
    FEl  = groupStates(FEl,vi,vpl,irow,icol,"fe"); % fast-equilibrium B.C. rates
    SSl  = groupStates(SSl,vi,vpl,irow,icol,"ss"); % steady state approx. rates
    TSTl = groupStates(TSTl,vi,vpl,irow,icol,"ts"); % TST rates
    NFl  = groupStates(NFl,vi,vpl,irow,icol,"nf"); % no-flux B.C. rates
    [irow,icol,know] = maxNZ(min(Kl,Kl'))
    vpl = vpln; vi = vin;
    nStates = size(Kl,1) % print progress
  endwhile
endfunction

function [Kln,vpl,vi] = groupStates(Kl,vi,vpl,is1,is2,kType)
% rate matrix reduction by 1 state
  is12 = min(is1,is2); is2 = max(is1,is2); is1 = is12;
  R12 = Kl(is1,is2); R21 = Kl(is2,is1);
  Kln = Kl; % new matrix is created, so that the gradual modifications do not affect
  Kl(is1,is2) = 0; Kl(is2,is1) = 0; Kl(is1,is1) = 0; Kl(is2,is2) = 0; % delete rates between s1 and s2, so they do not appear among neighbours
%  Kl(is1,find(Kl(is1,:)<-1e5))=0; Kl(is2,find(Kl(is2,:)<-1e5))=0; Kl(find(Kl(:,is1)<-1e5),is1)=0; Kl(find(Kl(:,is2)<-1e5),is1)=0; % discard too small rates
  vnei = unique([find(Kl(is1,:)),find(Kl(is2,:))]); % find all the neighbours of s1 and s2
  nj = size(vnei,2);
  if (nj>0)
    for j=vnei
      vpj = [vpl(is1);vpl(is2);vpl(j)]; vpj = vpj - logSumExp(vpj); % vector of 3 eq. popul.
      Kij = [0,R12,Kl(is1,j);R21,0,Kl(is2,j);Kl(j,is1),Kl(j,is2),0]; % 3x3 log rate matrix
      if (max(max(Kij))-min(min(Kij))>1e3)
	betweenObserv = [is1, is2, j]
	[kab,kba] = tstlog(Kij,vpj);
      elseif (prod(kType=="fe")==1)
	vKl1 = Kl(:,is1); vKl1(is2) = 0; vKl1(j) = 0; vKl1 = full(vKl1(find(vKl1))); kADl = logSumExp(vKl1); % adding the equilibration through the neighbours
	vKl2 = Kl(:,is2); vKl2(is1) = 0; vKl2(j) = 0; vKl2 = full(vKl2(find(vKl2))); kBDl = logSumExp(vKl2);
	sumflux = logSumExp([vpl(1)+kADl;vpl(2)+kBDl]); 
	f1l = vpl(1)+kADl-sumflux; 
	f2l = vpl(2)+kBDl-sumflux; 
	Kij(1,2) = logSumExp([Kij(1,2);kBDl+f1l]);
	Kij(2,1) = logSumExp([Kij(2,1);kADl+f2l]);
	[kab,kba] = analyt3x3(Kij,vpj);
      elseif (prod(kType=="nf")==1)
	[kab,kba] = analyt3x3(Kij,vpj);
      elseif (prod(kType=="ss")==1)
	[kab,kba] = sslog(Kij,vpj);
      else
	[kab,kba] = tstlog(Kij,vpj);
      endif
      Kln(is1,j) = kab;
      Kln(j,is1) = kba;
    endfor
  endif
  vpl(is1) = logSumExp([vpl(is1);vpl(is2)]); % grouped population
  vpl(is2) = []; Kln(is2,:) = []; Kln(:,is2) = []; % delete s2 from R and vp
  vi([find(vi==is1);find(vi==is2)]) = is1; % indices of all minima grouped into s1
  nStates = size(Kln,1);
  if (is2<=nStates)
    vi(find(vi>=is2))--; % deleting s2, R became smaller -> decrement all indices>=is2
  endif
endfunction
