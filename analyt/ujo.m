Kl = load('Kl').Kl;
vpl = load('vpl');
vil = load('vlump');
[vi,FEl,SSl,TSTl,NFl] = rrr(Kl,vpl,2,vil) % group into a 2-state model with arbitrary observables
