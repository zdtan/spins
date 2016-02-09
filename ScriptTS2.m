%% This creates thermal state in special way, not from random configuration
% starts with ground state for thin system (try constrainM = 0 versus constrainM = 1)
% copies the layers into full system
% heats up to low temperature
%% Includes the work of ScriptCorrelations

%% Really just for initialisation and running GenerateTS2
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%doMC = true;
%doSD = true;
%delta0 = 1;	% this is the value of delta for beta = 1
%numSweepMC = 10^4;
AFM = true;
sigma = 1;

% automatic initialisation
ID = ['a',num2str(alpha),'-b',num2str(beta),'-s',num2str(sigma),...
    '-',num2str(K),'x',num2str(L),'x',num2str(M),'_',UID];
% calculate adjacency matrix just once
Jmat = AdjacencyMatrix( alpha, K, L, M, AFM );

%% work on thin system
fakeM = 1;
JmatThin = AdjacencyMatrix( alpha, K, L, fakeM, AFM );
numSweepMCcool = 10^5;
betaCool = 10^3; 
numSweepSD = 10^5; 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Generate thermal state

	ticMC = tic;
	% random number generator setting for mixing things up!
	rng('shuffle')
	sinit = rng;
	[ spinMat, energy, acceptance ] = GenerateTS( ID, ...
		alpha, betaCool, sigma, K, L, fakeM, AFM, constrainM, JmatThin, sinit, numSweepMCcool, delta0 );
	timeMC = toc(ticMC);
	disp(['Time taken on Monte Carlo (thin) = ',num2str(timeMC./3600),' hrs']);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Generate GS using SD from thermal state, else load GS

	ticSD = tic;
	[ spinMat, energy ] = GenerateGS( ID, ...
		alpha, JmatThin, betaCool, sigma, K, L, fakeM, AFM, numSweepSD );
	timeSD = toc(ticSD);
	disp(['Time taken on Steepest Descent (thin) = ',num2str(timeSD./3600),' hrs']);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Generate GS for full system by copying
spinMat = repmat(spinMat,M,1); 

%% save spin configuration
spinfile = [ID,'/spins-a',num2str(alpha),'-SD', ...
    '-s',num2str(sigma), ...
    '-',num2str(K),'x',num2str(L),'x',num2str(M)];
if ~AFM
    spinfile = [spinfile,'_FM'];
end
save([spinfile,'.mat'],'spinMat',...
	'alpha','sigma','K','L','M','-mat');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% correlations for SD ground state
folder = [ID,'/correlationsSD'];
%% "create" directory for saving
[s,mess,messid] = mkdir(folder);
if s == 1
    % successful in creating directory
    if strcmp(messid,'MATLAB:MKDIR:DirectoryExists')
        % directory already exists
    else
        %disp(['== ABORT: directory does not previously exist! ==']);
        %return
    end
else
    disp(['== ABORT: failed to create directory! ==']);
    return
end

%% calculate correlations
CorrelationsAvg( spinMat, K, L, M, 'api-col', folder );

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Heat up to a thermal state

	ticMC2 = tic;
	% random number generator setting for mixing things up!
	rng('shuffle')
	sinit = rng;
	[ spinMat, energy, acc ] = GenerateTS2( spinMat, ID, ...
		alpha, beta, sigma, K, L, M, AFM, constrainM, Jmat, sinit, numSweepMC, delta0 );
	timeMC2 = toc(ticMC2);
	disp(['Time taken on Monte Carlo = ',num2str(timeMC./3600),' hrs']);

%	ticMC = tic;
	% random number generator setting for mixing things up!
%	rng('shuffle')
%	sinit = rng;

%delta = ;

%	[ spinMat, energy, acc ] = MonteCarlo( spinMat, Jmat, ...
%		alpha, beta, delta, sigma, K, L, M, numSweep, AFM, constrainM );

%	[ spinMat, energy, acceptance ] = GenerateTS( ID, ...
%		alpha, betaCool, sigma, K, L, fakeM, AFM, constrainM, JmatThin, sinit, numSweepMCcool, delta0 );
%	timeMC = toc(ticMC);
%	disp(['Time taken on Monte Carlo = ',num2str(timeMC./3600),' hrs']);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% correlations for thermal state
folder = [ID,'/correlations2'];
%% "create" directory for saving
[s,mess,messid] = mkdir(folder);
if s == 1
    % successful in creating directory
    if strcmp(messid,'MATLAB:MKDIR:DirectoryExists')
        % directory already exists
    else
        %disp(['== ABORT: directory does not previously exist! ==']);
        %return
    end
else
    disp(['== ABORT: failed to create directory! ==']);
    return
end

%% calculate correlations
CorrelationsAvg( spinMat, K, L, M, 'api-col', folder );

