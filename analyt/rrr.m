function [vi,FEl,SSl,TSTl,NFl] = rrr(Kl,vpl,nDesiredStates,vObserv,maxThresh)
% grouping minima stops when either a rate constant threshold or a number of states is reached
% returns matrices without diagonal terms
% grouping is determined by largest of the lowest rate constants calculated using no-flux BC
% Kl ... log of the rate matrix
% vpl ... log of equil populs
% nDesiredStates ... grouping ends when this number of states is reached
% vObserv ... determines the states that should not be grouped together, as they belong to different observables
  n = size(Kl,1); nStates = n; % the original number of states
  if (nargin<4) vObserv = ones(n,1); endif % if anything can be grouped, make vObserv vector
  if (nargin<5) maxThresh = -99e9; endif
  vObsInd = unique(vObserv); nObserv = size(vObsInd,1); % how many observable states there are
  if (nDesiredStates<nObserv) nDesiredStates = nObserv; endif % number of states at the end is at least the number of observables
  vi = [1:n]'; % index vector, always of size n
  rozdiel = max(max(Kl))+20
  Kl(find(Kl)) -= rozdiel;
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
  while and(nStates>nDesiredStates,know>maxThresh)
    [Kl,vpln,vin]   = groupStates(Kl,vi,vpl,irow,icol); % new decision matrix, assignment, and population vector
    [TSTl,FEl,NFl,SSl]  = groupStatesAll(TSTl,FEl,NFl,SSl,vi,vpl,irow,icol); % fast-equilibrium B.C. rates
    [irow,icol,know] = maxNZ(min(Kl,Kl'))
    kTSl = TSTl(irow,icol)+rozdiel
    kFEl = FEl(irow,icol)+rozdiel
    kNFl = NFl(irow,icol)+rozdiel
    kSSl = SSl(irow,icol)+rozdiel
    vpl = vpln; vi = vin;
    nStates = size(Kl,1) % print progress
  endwhile
  FEl += rozdiel;
  NFl += rozdiel;
  TSTl += rozdiel;
  SSl += rozdiel;
endfunction

function [TSTl,FEl,NFl,SSl]  = groupStatesAll(TSTl,FEl,NFl,SSl,vi,vpl,is1,is2)
% rate matrix reduction by 1 state
% all types of rate definitions
  nic = -9e9;
%  TSTln = TSTl; FEln = FEl; NFln = NFl; SSln = SSl;
% is1 will be the index after the rearrangement
  is12 = min(is1,is2); is2 = max(is1,is2); is1 = is12;
% for FEl all will change
  vFElis1j = FEl(is1,:); vFEljis1 = FEl(:,is1);
  vFElis2j = FEl(is2,:); vFEljis2 = FEl(:,is2);
% first find neighbours - they must be the same for all matrices -> use TST matrix
% neighbours of is1 except is2
  vil1 = setdiff(intersect(find(TSTl(is1,:)>-1e8),find(TSTl(is1,:))),[is1,is2]);
%  vTSTis1jrates = TSTl(is1,vil1)
%  vFEis1jrates = vFElis1j([is2,vil1])
%  vFEjis1rates = vFEljis1([is2,vil1])
% neighbours of is2 except is1
  vil2 = setdiff(intersect(find(TSTl(is2,:)>-1e8),find(TSTl(is2,:))),[is1,is2]);
%  vFEis2jrates = vFElis2j([is1,vil2])
%  vFEjis2rates = vFEljis2([is1,vil2])
  vnei = unique([vil1,vil2]);
  [nnj,nj] = size(vnei); % if both are 
  if or(nj==0,nnj==0)
    there_were_no_neighbours_wtf = [is1, is2]
  else
    for j=vnei
      vpj = [vpl(is1);vpl(is2);vpl(j)]; vpj = vpj - logSumExp(vpj); % vector of 3 eq. popul.; normalis can make something 0
% do TST first
%      ktyp = "tst"
      TSTij = [nic,           TSTl(is1,is2), TSTl(is1,j);
	       TSTl(is2,is1),           nic, TSTl(is2,j);
	       TSTl(  j,is1), TSTl(  j,is2),        nic];
      [tstlkab,tstlkba] = tstlog(TSTij,vpj)
      TSTl(is1,j) = tstlkab; TSTl(j,is1) = tstlkba;
% NF
%      ktyp = "nf"
      NFij = [nic,          NFl(is1,is2), NFl(is1,j);
	      NFl(is2,is1),          nic, NFl(is2,j);
	      NFl(  j,is1), NFl(  j,is2),       nic];
      [nflkab,nflkba] = analyt3x3(NFij,vpj)
      NFl(is1,j) = nflkab; NFl(j,is1) = nflkba;
% FE
%      ktyp = "fe"
      FEij = [nic,           vFElis1j(is2), vFElis1j(j);
	      vFElis2j(is1),           nic, vFElis2j(j);
	        vFEljis1(j),   vFEljis2(j),       nic];
      vil1j = setdiff(vil1,j); [nun1,nun2] = size(vil1j); 
      if (nun1*nun2>0) kDAl = logSumExp(vFEljis1(vil1j)) else kDAl = nic; endif
      vil2j = setdiff(vil2,j); [nun1,nun2] = size(vil2j); 
      if (nun1*nun2>0) kDBl = logSumExp(vFEljis2(vil2j)) else kDBl = nic; endif
      sumflux = logSumExp([vpj(1)+kDAl;vpj(2)+kDBl]);
      f1l = vpj(1)+kDAl-sumflux;
      f2l = vpj(2)+kDBl-sumflux;
      FEij(1,2) = logSumExp([FEij(1,2);kDBl+f1l]);
      FEij(2,1) = logSumExp([FEij(2,1);kDAl+f2l]);
      [felkab,felkba] = analyt3x3(FEij,vpj)
      FEl(is1,j) = felkab; FEl(j,is1) = felkba;
% SS
%      ktyp = "ss"
      SSij = [nic,          SSl(is1,is2), SSl(is1,j);
	      SSl(is2,is1),          nic, SSl(is2,j);
	      SSl(  j,is1), SSl(  j,is2),       nic];
      [sslkab,sslkba] = sslog(SSij,vpj)
      SSl(is1,j) = sslkab; SSl(j,is1) = sslkba;
    endfor
    Porovn = [TSTl(is1,vnei)',FEl(is1,vnei)',NFl(is1,vnei)',SSl(is1,vnei)']
  endif
  TSTl(is2,:) = []; TSTl(:,is2) = [];
  FEl(is2,:) = []; FEl(:,is2) = [];
  NFl(is2,:) = []; NFl(:,is2) = [];
  SSl(is2,:) = []; SSl(:,is2) = [];
endfunction

function [Kln,vpl,vi] = groupStates(Kl,vi,vpl,is1,is2)
% rate matrix reduction by 1 state
% the decision matrix only
  malo = -9e9;
  is12 = min(is1,is2); is2 = max(is1,is2); is1 = is12;
  R12 = Kl(is1,is2); R21 = Kl(is2,is1);
  Kln = Kl; % new matrix is created, so that the gradual modifications do not affect
  Kl(is1,is2) = malo; Kl(is2,is1) = malo; Kl(is1,is1) = malo; Kl(is2,is2) = malo;
  vnei = unique([intersect(find(Kl(is1,:)>-1e8),find(Kl(is1,:))),intersect(find(Kl(is2,:)>-1e8),find(Kl(is2,:)))]);
  [nnj,nj] = size(vnei)
  if and(nj>0,nnj>0)
    for j=vnei
      vpj = [vpl(is1);vpl(is2);vpl(j)]; vpj = vpj - logSumExp(vpj); % vector of 3 eq. popul.
      Kij = [malo,R12,Kl(is1,j);R21,malo,Kl(is2,j);Kl(j,is1),Kl(j,is2),malo]; % 3x3 log rate matrix
      [kab,kba] = analyt3x3(Kij,vpj);
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

