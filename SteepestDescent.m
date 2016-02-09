%% taken from 20140826-spinwaves

function [ spinMat, energy ] = SteepestDescent( initialSpins, Jmat, ...
    alpha, delta, sigma, K, L, M, numSweep, AFM )
%SteepestDescent - Goes to every site and aligns spin exactly with local field for
%       delta = 1, or under/over-shoots for delta less than or more than 1.
% all sites are visited in repeatably random fashion.
% [this is done elsewhere] Save in folder 'ID' (which is already present)
%%%%%%%%%%%%%%%%%%%%%%%%%

% setup random number generator
rng('default');	% equivalent to rng(0,'twister');

%%% create directory for saving
%[s,mess,messid] = mkdir(ID);
%if s == 1
%    % successful in creating directory
%    if strcmp(messid,'MATLAB:MKDIR:DirectoryExists')
%        %~ disp(['== WARNING: directory already exists! ==']);
%        %~ return
%    else
%        disp(['== WARNING: directory not created yet! ==']);
%        return
%    end
%else
%    disp(['== ABORT: failed to create directory! ==']);
%    return
%end

%% system size
N = 8*K*L*M; % 8 sites per unit cell

%% initialise result
spinMat = initialSpins;

%% system energy per spin
energy = zeros(1,numSweep+1);

% calculate energy of system at the very start
energy(1) = 0.5*sum(dot(spinMat,Jmat*spinMat,2))/N ...
	- (-3/4 - (alpha*sigma)^2/2); %...
	%+ constrainM*sum(abs(totalM).^2)/N;

%% spinVec modification in Steepest Descent move.
% normalise fieldX to the length S of spin1, call it spinX
% spin2 = spin1 + delta*(spinX-spin1)
%  then scale back final answer to length S
function spin2 = shiftSpin(spin1,delta,fieldX) % fieldX is reference field
	spinLength = norm(spin1,2);
	spinX = fieldX.*spinLength./norm(fieldX,2);
	spin2 = spin1 + delta*(spinX-spin1);
	spin2 = spin2.*spinLength./norm(spin2,2);
end

%% begin sweeps
disp(['SET alpha = ',num2str(alpha),' , sigma = ',num2str(sigma), ...
    ', for ', ...
    num2str(numSweep),' steepest descent runs...']);
for t = 1:numSweep
    % disp(['...Steepest Descent run ',num2str(t),' of ',num2str(numSweep),'...']);
    % tic;
    energy(t+1) = energy(t); % initialise energy from previous sweep

    %% visit each site once
    visitOrder = randperm(N);	% visit all sites in random fashion
    for s = 1:N
        % find spin, find local field (of same length as spin), then move.
        n = visitOrder(s);
        spin1 = spinMat(n,:);
        localField = -Jmat(n,:)*spinMat;
        spin2 = shiftSpin(spin1,delta,localField);
        % turn spin and update energy
        spinMat(n,:) = spin2;
        energy(t+1) = energy(t+1) - dot((spin2-spin1),localField,2)/N;
    end

%     %% check energy calculation
%	 energycheck(t+1) = 0.5*sum(dot(spinMat,Jmat*spinMat,2))/N ...
%		- (-3/4 - (alpha*sigma)^2/2);
%     % [old version] calculate total energy (account for double counting)
%     for n = 1:N
%         energycheck(t+1) = energycheck(t+1) ...
%             - 0.5*(spins(n)*LocalField(n,K,L,M,spins,alpha,AFM));
%     end
end

end
