%% Really just for initialisation and running GenerateTS
% does simulated annealing if doMC is true
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%delta0 = 1;	% this is the value of delta for beta = 1
doMC = true;
%numSweepMC = 10^4;
AFM = true;
sigma = 1;

% automatic initialisation
ID = ['a',num2str(alpha),'-b',num2str(beta),'-s',num2str(sigma),...
    '-',num2str(K),'x',num2str(L),'x',num2str(M),'_',UID];
% calculate adjacency matrix just once
Jmat = AdjacencyMatrix( alpha, K, L, M, AFM );

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Generate thermal state
if doMC
	ticMC = tic;
	% random number generator setting for mixing things up!
	rng('shuffle')
	sinit = rng;
	[ spinMat, energy, acceptance ] = GenerateTS( ID, ...
		alpha, beta, sigma, K, L, M, AFM, constrainM, Jmat, sinit, numSweepMC, delta0 );
	timeMC = toc(ticMC);
	disp(['Time taken on Monte Carlo = ',num2str(timeMC./3600),' hrs']);
end

